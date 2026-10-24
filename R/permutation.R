#' Permutation-based network proximity (for validation and comparison)
#'
#' Degree-bin-preserving Monte-Carlo permutation version of the proximity
#' test. Used to validate the exact null distributions of
#' \code{\link{netEnrich}} and to compute the symmetric module-to-module
#' distance used by NetFunrankR v1 (Menche et al. 2015) or the closest
#' distance of Guney et al. (2016) with a z-score.
#'
#' @inheritParams netEnrich
#' @param nPerm number of permutations.
#' @param statistic \code{"closest"} (\eqn{d_c}) or \code{"symmetric"}
#'   (\eqn{d_{XY}}).
#' @param batch permutations processed together in one multi-source search
#'   (symmetric statistic only).
#' @returns \code{DataFrame} with observed statistic, null mean and
#'   standard deviation, z-score, empirical p-value
#'   \eqn{(1 + \#\{null \le obs\}) / (1 + B)} and normal-approximation
#'   p-value.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' kegg <- loadPathways("KEGG")[1:20]
#' set.seed(1)
#' head(netEnrichPerm(COPD_proteinCodingGenes, kegg, net, nPerm = 200))
netEnrichPerm <- function(query, geneSets, network, distances = NULL,
    nPerm = 1000L, statistic = c("closest", "symmetric"), minSize = 5L,
    maxSize = Inf, batch = 50L) {
    statistic <- match.arg(statistic)
    geneSets <- asGeneSetList(geneSets)
    if (is.null(distances))
        distances <- geneSetDistances(network, geneSets)
    nIn <- colSums(distances == 0L, na.rm = TRUE)
    D <- distances[, nIn >= minSize & nIn <= maxSize, drop = FALSE]
    qi <- which(nodeNames(network) %in% query)
    draw <- .binnedSampler(qi, network)
    null <- if (statistic == "closest")
        vapply(seq_len(nPerm), function(i) colMeans(D[draw(), , drop = FALSE]),
            numeric(ncol(D)))
    else .symmetricNull(D, network, draw, nPerm, batch)
    obs <- if (statistic == "closest") colMeans(D[qi, , drop = FALSE]) else
        .symmetricStat(D, qi, geneSetDistances(network,
            list(q = nodeNames(network)[qi]))[, 1])
    null <- matrix(null, nrow = ncol(D))
    mu <- rowMeans(null)
    sdv <- apply(null, 1, stats::sd)
    z <- (obs - mu) / sdv
    S4Vectors::DataFrame(set = colnames(D), statistic = statistic,
        observed = obs, null_mean = mu, null_sd = sdv, z = z,
        p_perm = (1 + rowSums(null <= obs)) / (1 + nPerm),
        p_norm = stats::pnorm(z))
}

## returns a function drawing degree-bin matched random node indices
.binnedSampler <- function(qi, network) {
    bin <- degreeBin(network)
    pools <- split(seq_along(bin), bin)
    qbin <- bin[qi]
    function() {
        r <- integer(length(qi))
        for (b in unique(qbin)) {
            w <- which(qbin == b)
            pool <- pools[[as.character(b)]]
            r[w] <- pool[sample.int(length(pool), length(w))]
        }
        r
    }
}

.symmetricStat <- function(D, xi, dX) {
    member <- D == 0L
    (colSums(D[xi, , drop = FALSE]) + colSums(member * dX)) /
        (length(xi) + colSums(member))
}

.symmetricNull <- function(D, network, draw, nPerm, batch) {
    genes <- nodeNames(network)
    null <- matrix(NA_real_, ncol(D), nPerm)
    for (st in seq(1L, nPerm, by = batch)) {
        ids <- seq.int(st, min(nPerm, st + batch - 1L))
        R <- lapply(ids, function(i) draw())
        dR <- geneSetDistances(network,
            stats::setNames(lapply(R, function(r) genes[r]), ids))
        for (j in seq_along(ids))
            null[, ids[j]] <- .symmetricStat(D, R[[j]], dR[, j])
    }
    null
}

#' Over-representation analysis (hypergeometric test)
#'
#' One-sided test \eqn{P(X \ge x)} of the overlap between query genes and
#' each gene set against a common background, as in DAVID or
#' clusterProfiler. Provided as the reference method.
#'
#' @inheritParams netEnrich
#' @param universe background genes; default: all genes annotated in
#'   \code{geneSets}.
#' @returns \code{DataFrame} ordered by p-value.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' head(oraTest(COPD_proteinCodingGenes, loadPathways("KEGG")))
oraTest <- function(query, geneSets, universe = NULL, minSize = 5L,
    maxSize = Inf, pAdjust = "BH") {
    geneSets <- asGeneSetList(geneSets)
    if (is.null(universe))
        universe <- unique(unlist(geneSets, use.names = FALSE))
    q <- intersect(unique(query), universe)
    inU <- lapply(geneSets, intersect, universe)
    ov <- lapply(inU, intersect, q)
    M <- unname(lengths(inU))
    x <- unname(lengths(ov))
    res <- S4Vectors::DataFrame(set = names(geneSets), set_size = M,
        overlap = x, expected = length(q) * M / length(universe),
        pvalue = stats::phyper(x - 1L, M, length(universe) - M, length(q),
            lower.tail = FALSE),
        overlap_genes = unname(vapply(ov, paste, character(1),
            collapse = "/")))
    res <- res[M >= minSize & M <= maxSize, ]
    res$padj <- stats::p.adjust(res$pvalue, pAdjust)
    res <- res[order(res$pvalue), ]
    rownames(res) <- NULL
    res
}

#' Test many query lists against one library
#'
#' Computes the gene-set distances once and runs \code{\link{netEnrich}}
#' for every query list, optionally in parallel.
#'
#' @inheritParams netEnrich
#' @param queries named list of character vectors.
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}.
#' @param ... further arguments to \code{\link{netEnrich}}.
#' @returns named list of \code{\linkS4class{NFResult}} objects.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' kegg <- loadPathways("KEGG")
#' q <- list(full = COPD_proteinCodingGenes,
#'     half = COPD_proteinCodingGenes[1:150])
#' res <- netEnrichBatch(q, kegg, net)
#' vapply(res, function(r) sum(enrichResults(r)$fdr_combined < 0.05),
#'     integer(1))
netEnrichBatch <- function(queries, geneSets, network, ...,
    BPPARAM = BiocParallel::SerialParam()) {
    if (!is.list(queries) || is.null(names(queries)))
        stop("'queries' must be a named list of gene vectors")
    geneSets <- asGeneSetList(geneSets)
    D <- geneSetDistances(network, geneSets)
    BiocParallel::bplapply(queries, netEnrich, geneSets = geneSets,
        network = network, distances = D, ..., BPPARAM = BPPARAM)
}

#' Exact null distributions for one gene set
#'
#' Returns the null probability mass functions of the overlap, total
#' distance and beyond-overlap distance used by \code{\link{netEnrich}}
#' for a single gene set (exact convolution), e.g. for plotting with
#' \code{\link{plotNullDistribution}}.
#'
#' @inheritParams netEnrich
#' @param set name of the gene set.
#' @returns list with the observed values (\code{n}, \code{overlap},
#'   \code{S_obs}, \code{d_query}) and pmfs \code{O}, \code{S}, \code{Sb}.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' kegg <- loadPathways("KEGG")
#' nd <- nullDistribution(COPD_proteinCodingGenes, kegg, net,
#'     "Nicotine addiction")
#' nd$S_obs / nd$n
nullDistribution <- function(query, geneSets, network, set) {
    geneSets <- asGeneSetList(geneSets)[set]
    if (is.null(geneSets[[1]]))
        stop("gene set '", set, "' not found")
    d <- geneSetDistances(network, geneSets)[, 1]
    prep <- .prepareQuery(query, network, NULL)
    kern <- .binKernels(d, prep)
    member <- d[prep$qi] == 0L
    pm <- exact_null_cpp(kern$f, prep$mB,
        tabulate(prep$qbin[!member], nbins = prep$B))
    c(pm, list(n = length(prep$qi), overlap = sum(member),
        S_obs = sum(d[prep$qi]),
        d_query = stats::setNames(d[prep$qi], prep$genes[prep$qi])))
}

#' Decomposed network-proximity enrichment test
#'
#' Tests whether query genes are closer to each gene set in the interactome
#' than degree-matched random genes, and decomposes the evidence into
#' (i) direct membership (overlap) and (ii) interactome proximity of the
#' non-member query genes ("beyond overlap"). All null distributions are
#' exact (discrete convolution over degree bins; saddlepoint approximation
#' for long lists), so no permutations are needed.
#'
#' Under the null hypothesis each query gene is replaced by a random gene of
#' the same degree bin (drawn from \code{background} if supplied). The
#' following one-sided tests are reported for every gene set \eqn{Y}:
#' \describe{
#'   \item{\code{p_overlap}}{degree-aware overlap test, \eqn{P(O \ge o)} for
#'     a Poisson-binomial overlap \eqn{O}.}
#'   \item{\code{p_linker}}{number of non-member query genes adjacent to
#'     \eqn{Y} (linker genes), conditional on the membership pattern.}
#'   \item{\code{p_beyond}}{summed distance of non-member query genes to
#'     \eqn{Y}, conditional on the membership pattern.}
#'   \item{\code{p_proximity}}{total closest-distance proximity
#'     \eqn{d_c} (Guney et al. 2016).}
#'   \item{\code{p_combined}}{Fisher combination of the mid-p values of
#'     \code{p_overlap} and the beyond-overlap test chosen by
#'     \code{beyond}; valid because the two are independent under the
#'     null.}
#'   \item{\code{p_hyper}}{classical hypergeometric ORA, for reference.}
#' }
#' Each p-value is adjusted across gene sets (\code{fdr_*} columns).
#'
#' @param query character vector of query genes (same identifiers as the
#'   network).
#' @param geneSets gene sets in any format accepted by
#'   \code{\link{asGeneSetList}}.
#' @param network an \code{\linkS4class{NFNetwork}}.
#' @param distances optional result of \code{\link{geneSetDistances}} for
#'   \code{network} and \code{geneSets}; supply it when testing several
#'   query lists against the same library.
#' @param minSize,maxSize gene sets with fewer/more members in the network
#'   are not tested.
#' @param universe background for the hypergeometric ORA column; default:
#'   all genes annotated in \code{geneSets}.
#' @param background optional genes that could have been selected (e.g. the
#'   genes measured in an expression experiment); null genes are then drawn
#'   from these only.
#' @param beyond beyond-overlap test entering \code{p_combined}:
#'   \code{"linker"} (default) or \code{"distance"}.
#' @param exactMax lists with at most this many network genes use exact
#'   convolution, longer lists a saddlepoint approximation.
#' @param pAdjust multiple-testing method passed to
#'   \code{\link[stats]{p.adjust}}.
#' @returns an \code{\linkS4class{NFResult}}; rows of
#'   \code{enrichResults()} are ordered by \code{p_combined}.
#' @references Guney E, et al. (2016) Network-based in silico drug efficacy
#'   screening. \emph{Nature Communications} 7:10331.
#' @seealso \code{\link{evidenceClass}}, \code{\link{plotEvidenceMap}},
#'   \code{\link{netEnrichBatch}}
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' kegg <- loadPathways("KEGG")
#' res <- netEnrich(COPD_proteinCodingGenes, kegg, net)
#' res
#' head(enrichResults(res)[, c("set", "overlap", "n_linker", "p_combined")])
netEnrich <- function(query, geneSets, network, distances = NULL,
    minSize = 5L, maxSize = Inf, universe = NULL, background = NULL,
    beyond = c("linker", "distance"), exactMax = 500L, pAdjust = "BH") {
    beyond <- match.arg(beyond)
    stopifnot(methods::is(network, "NFNetwork"))
    geneSets <- asGeneSetList(geneSets)
    query <- unique(as.character(query))
    if (is.null(distances))
        distances <- geneSetDistances(network, geneSets)
    if (!identical(rownames(distances), nodeNames(network)))
        stop("'distances' was not computed for this network")
    geneSets <- geneSets[colnames(distances)]
    prep <- .prepareQuery(query, network, background)
    if (is.null(universe))
        universe <- unique(unlist(geneSets, use.names = FALSE))
    nIn <- colSums(distances == 0L, na.rm = TRUE)
    keep <- which(nIn >= minSize & nIn <= maxSize)
    if (!length(keep))
        stop("no gene set has between minSize and maxSize network genes")
    rows <- lapply(keep, function(k) {
        .testOneSet(distances[, k], prep, beyond, exactMax,
            geneSets[[k]], universe, query)
    })
    res <- S4Vectors::DataFrame(do.call(rbind, rows))
    res$set <- colnames(distances)[keep]
    res$set_size <- lengths(geneSets[keep])
    res$set_in_net <- unname(nIn[keep])
    res <- .finaliseResults(res, pAdjust)
    methods::new("NFResult", results = res,
        query = nodeNames(network)[prep$qi],
        unmatched = setdiff(query, nodeNames(network)),
        geneSets = geneSets[res$set],
        params = list(beyond = beyond, minSize = minSize, maxSize = maxSize,
            exactMax = exactMax, pAdjust = pAdjust,
            network = network@source, background = !is.null(background)))
}

.prepareQuery <- function(query, network, background) {
    genes <- nodeNames(network)
    qi <- which(genes %in% query)
    if (length(qi) < 2L)
        stop("fewer than 2 query genes are in the network")
    bin <- degreeBin(network)
    B <- max(bin)
    bg <- if (is.null(background)) seq_along(genes) else
        sort(union(which(genes %in% background), qi))
    list(qi = qi, qbin = bin[qi], B = B, mB = tabulate(bin[qi], nbins = B),
        bg = bg, bgBin = bin[bg], genes = genes)
}

## Bin-specific kernels: f (distance), fO (membership), g (distance of a
## non-member), fL (linker indicator of a non-member).
.binKernels <- function(d, prep) {
    B <- prep$B
    dmax <- max(d, na.rm = TRUE)
    cnt <- matrix(tabulate(d[prep$bg] * B + prep$bgBin,
        nbins = B * (dmax + 1L)), nrow = B, ncol = dmax + 1L)
    f <- cnt / pmax(rowSums(cnt), 1)
    g <- f
    g[, 1L] <- 0
    rs <- rowSums(g)
    g <- g / pmax(rs, 1e-300)
    g[rs == 0, 2L] <- 1
    list(f = f, fO = cbind(1 - f[, 1L], f[, 1L]), g = g,
        fL = cbind(1 - g[, 2L], g[, 2L]))
}

.testOneSet <- function(d, prep, beyond, exactMax, setGenes, universe,
    query) {
    kern <- .binKernels(d, prep)
    dq <- d[prep$qi]
    member <- dq == 0L
    obs <- list(o = sum(member), s = sum(dq), l = sum(dq == 1L),
        mNon = tabulate(prep$qbin[!member], nbins = prep$B),
        anyNon = any(!member))
    p <- if (length(prep$qi) <= exactMax) .exactTails(kern, prep, obs) else
        .saddleTails(kern, prep, obs)
    pc <- if (beyond == "linker") p$linker else p$beyond
    chi <- -2 * (log(max(p$overlap[["mid"]], 1e-300)) +
        log(max(pc[["mid"]], 1e-300)))
    momS <- .sumMoments(kern$f, prep$mB)
    momB <- .sumMoments(kern$g, obs$mNon)
    nq <- length(prep$qi)
    nNon <- nq - obs$o
    data.frame(
        n_query_net = nq, overlap = obs$o,
        overlap_expected = sum(prep$mB * kern$f[, 1L]),
        n_linker = obs$l,
        linker_expected = sum(obs$mNon * kern$g[, 2L]),
        d_c = obs$s / nq,
        d_c_expected = momS[["mean"]] / nq,
        z_proximity = (momS[["mean"]] - obs$s) / max(momS[["sd"]], 1e-12),
        z_beyond = if (nNon) (momB[["mean"]] - obs$s) /
            max(momB[["sd"]], 1e-12) else NA_real_,
        p_hyper = .hyperP(query, setGenes, universe),
        p_overlap = p$overlap[["p"]], p_linker = p$linker[["p"]],
        p_beyond = p$beyond[["p"]], p_proximity = p$proximity[["p"]],
        p_combined = stats::pchisq(chi, df = 4, lower.tail = FALSE),
        overlap_genes = paste(prep$genes[prep$qi][member], collapse = "/"),
        linker_genes = paste(prep$genes[prep$qi][dq == 1L], collapse = "/"))
}

.exactTails <- function(kern, prep, obs) {
    nl <- exact_null_cpp(kern$f, prep$mB, obs$mNon)
    one <- c(p = 1, mid = 1)
    list(overlap = .upperP(nl$O, obs$o),
        proximity = .lowerP(nl$S, obs$s),
        beyond = if (obs$anyNon) .lowerP(nl$Sb, obs$s) else one,
        linker = if (obs$anyNon) .upperP(exact_null_cpp(
            kern$fL[, c(2L, 1L), drop = FALSE], obs$mNon,
            integer(prep$B))$O, obs$l) else one)
}

.saddleTails <- function(kern, prep, obs) {
    up <- function(f, m, k) {
        a <- .spUpper(f, m, k)
        c(p = a, mid = (a + .spUpper(f, m, k + 1L)) / 2)
    }
    lo <- function(f, m, k) {
        a <- .spLower(f, m, k)
        c(p = a, mid = (a + .spLower(f, m, k - 1L)) / 2)
    }
    one <- c(p = 1, mid = 1)
    list(overlap = up(kern$fO, prep$mB, obs$o),
        proximity = lo(kern$f, prep$mB, obs$s),
        beyond = if (obs$anyNon) lo(kern$g, obs$mNon, obs$s) else one,
        linker = if (obs$anyNon) up(kern$fL, obs$mNon, obs$l) else one)
}

## Tail probabilities of a pmf on 0..(length - 1) with mid-p values
.lowerP <- function(p, s) {
    cp <- cumsum(p)
    i <- s + 1L
    c(p = min(1, cp[i]), mid = min(1, cp[i] - p[i] / 2))
}
.upperP <- function(p, o) {
    tail <- sum(p[(o + 1L):length(p)])
    c(p = min(1, tail), mid = min(1, tail - p[o + 1L] / 2))
}

.hyperP <- function(query, setGenes, universe) {
    q <- intersect(query, universe)
    M <- length(intersect(setGenes, universe))
    x <- length(intersect(q, setGenes))
    stats::phyper(x - 1L, M, length(universe) - M, length(q),
        lower.tail = FALSE)
}

.finaliseResults <- function(res, pAdjust) {
    first <- c("set", "set_size", "set_in_net")
    res <- res[, c(first, setdiff(colnames(res), first))]
    for (p in c("p_hyper", "p_overlap", "p_linker", "p_beyond",
        "p_proximity", "p_combined"))
        res[[sub("^p_", "fdr_", p)]] <- stats::p.adjust(res[[p]], pAdjust)
    res <- res[order(res$p_combined, res$p_proximity), ]
    rownames(res) <- NULL
    res
}

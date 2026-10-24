#' Standardise gene sets to a named list
#'
#' Accepts the common gene-set representations used in R and Bioconductor.
#'
#' @param geneSets one of: a named list of character vectors; a two-column
#'   \code{data.frame} (gene, set), e.g. the bundled Enrichr libraries with
#'   columns \code{GeneSymbol} and \code{PathwayName}; a
#'   \code{GSEABase::GeneSetCollection}; or the path to a GMT file.
#' @returns named list of unique character vectors.
#' @export
#' @examples
#' data(KEGG_2021_Human, package = "NetFunrankR")
#' asGeneSetList(list(A = c("TP53", "MDM2"), B = "EGFR"))
#' head(names(asGeneSetList(KEGG_2021_Human)))
asGeneSetList <- function(geneSets) {
    if (is.character(geneSets) && length(geneSets) == 1L &&
        file.exists(geneSets))
        return(readGMT(geneSets))
    if (methods::is(geneSets, "GeneSetCollection")) {
        if (!requireNamespace("GSEABase", quietly = TRUE))
            stop("install 'GSEABase' to use GeneSetCollection objects")
        geneSets <- GSEABase::geneIds(geneSets)
    }
    if (is.data.frame(geneSets))
        geneSets <- .tableToList(geneSets)
    if (!is.list(geneSets) || is.null(names(geneSets)) ||
        anyNA(names(geneSets)))
        stop("'geneSets' must be a named list, a two-column data.frame, ",
            "a GeneSetCollection or a GMT file")
    lapply(geneSets, function(x) unique(as.character(x[!is.na(x)])))
}

.tableToList <- function(df) {
    gcol <- intersect(c("GeneSymbol", "gene", "Gene", "gene_symbol"),
        colnames(df))[1]
    scol <- intersect(c("PathwayName", "term", "Term", "set", "gs_name"),
        colnames(df))[1]
    if (is.na(gcol) || is.na(scol)) {
        gcol <- colnames(df)[1]
        scol <- colnames(df)[2]
    }
    split(as.character(df[[gcol]]), as.character(df[[scol]]))
}

#' Read gene sets from a GMT file
#'
#' @param file path to a GMT file (one set per line: name, description,
#'   genes; tab-separated).
#' @returns named list of character vectors.
#' @export
#' @examples
#' gmt <- tempfile(fileext = ".gmt")
#' writeLines(c("SET1\tdesc\tTP53\tMDM2", "SET2\tdesc\tEGFR\tERBB2"), gmt)
#' readGMT(gmt)
readGMT <- function(file) {
    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(lines)]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    sets <- lapply(parts, function(p) unique(p[-c(1, 2)][nzchar(p[-c(1, 2)])]))
    names(sets) <- vapply(parts, `[`, character(1), 1)
    sets
}

#' Load a bundled pathway library
#'
#' @param db one of \code{"KEGG"}, \code{"Reactome"} or
#'   \code{"WikiPathways"} (Enrichr libraries, HGNC symbols).
#' @returns named list of gene-symbol vectors.
#' @seealso \code{\link{KEGG_2021_Human}}
#' @export
#' @examples
#' kegg <- loadPathways("KEGG")
#' length(kegg)
loadPathways <- function(db = c("KEGG", "Reactome", "WikiPathways")) {
    db <- match.arg(db)
    obj <- switch(db,
        KEGG = "KEGG_2021_Human",
        Reactome = "Reactome_Pathways_2024",
        WikiPathways = "WikiPathways_2024_Human")
    env <- new.env(parent = emptyenv())
    utils::data(list = obj, package = "NetFunrankR", envir = env)
    asGeneSetList(get(obj, envir = env))
}

#' Shortest-path distance from every node to every gene set
#'
#' Computes \eqn{d_Y(v) = \min_{y \in Y} d(v, y)} for all nodes \eqn{v}
#' and gene sets \eqn{Y} at once by a multi-source breadth-first search
#' expressed as one sparse matrix product per search layer. Compute it once
#' per network and library and pass it to \code{\link{netEnrich}} when
#' testing several query lists.
#'
#' @param network an \code{\linkS4class{NFNetwork}}.
#' @param geneSets gene sets in any format accepted by
#'   \code{\link{asGeneSetList}}.
#' @param chunkSize number of gene sets searched together; smaller values
#'   lower the peak memory use (useful on small servers), larger values are
#'   slightly faster.
#' @returns integer matrix (nodes x gene sets); \code{NA} for sets with no
#'   gene in the network.
#' @export
#' @examples
#' net <- stringNetwork()
#' D <- geneSetDistances(net, loadPathways("KEGG")[1:5])
#' table(D[, 1])
geneSetDistances <- function(network, geneSets, chunkSize = 200L) {
    geneSets <- asGeneSetList(geneSets)
    genes <- nodeNames(network)
    idx <- lapply(geneSets, function(s) which(genes %in% s))
    D <- matrix(NA_integer_, length(genes), length(idx),
        dimnames = list(genes, names(geneSets)))
    adj <- adjacency(network)
    blocks <- split(seq_along(idx), ceiling(seq_along(idx) / chunkSize))
    for (b in blocks)
        D[, b] <- .bfsBlock(adj, idx[b], length(genes))
    D
}

## Multi-source breadth-first search for a block of gene sets.
.bfsBlock <- function(adj, idx, n) {
    D <- matrix(NA_integer_, n, length(idx))
    frontier <- matrix(FALSE, n, length(idx))
    for (k in seq_along(idx)) frontier[idx[[k]], k] <- TRUE
    D[frontier] <- 0L
    visited <- frontier
    step <- 0L
    while (any(frontier)) {
        step <- step + 1L
        frontier <- as.matrix(adj %*% (frontier * 1)) > 0 & !visited
        D[frontier] <- step
        visited <- visited | frontier
    }
    D
}

#' Build an interactome for network-proximity enrichment
#'
#' Converts an undirected edge list into an \code{\linkS4class{NFNetwork}}:
#' the largest connected component (LCC) of the simple graph, its sparse
#' adjacency matrix, node degrees and degree bins.
#'
#' Degree bins follow Guney et al. (2016): nodes are sorted by degree and
#' consecutive degree values are pooled until a bin holds at least
#' \code{minBinSize} nodes. Under the null model, each query gene is
#' replaced by a random gene of its own bin.
#'
#' @param edges a \code{data.frame}, matrix or \code{igraph} object; for
#'   tables the first two columns are the interacting node IDs.
#' @param lcc logical; restrict to the largest connected component
#'   (default \code{TRUE}, needed for finite shortest-path distances).
#' @param minBinSize minimum number of nodes per degree bin.
#' @param source short description of the network source.
#' @returns an \code{\linkS4class{NFNetwork}}.
#' @seealso \code{\link{stringNetwork}} for the bundled STRING interactome.
#' @references Guney E, Menche J, Vidal M, Barabasi AL (2016).
#'   Network-based in silico drug efficacy screening.
#'   \emph{Nature Communications} 7:10331.
#' @export
#' @examples
#' edges <- data.frame(from = c("A", "B", "C", "D"),
#'     to = c("B", "C", "D", "A"))
#' buildNetwork(edges, minBinSize = 2)
buildNetwork <- function(edges, lcc = TRUE, minBinSize = 100L,
    source = "user-supplied") {
    if (inherits(edges, "igraph"))
        edges <- igraph::as_data_frame(edges, what = "edges")
    edges <- as.data.frame(edges)[, c(1, 2)]
    colnames(edges) <- c("from", "to")
    edges$from <- as.character(edges$from)
    edges$to <- as.character(edges$to)
    ok <- !is.na(edges$from) & !is.na(edges$to) & edges$from != edges$to
    if (!any(ok))
        stop("no valid edges supplied")
    g <- igraph::simplify(
        igraph::graph_from_data_frame(edges[ok, ], directed = FALSE))
    if (lcc) {
        comp <- igraph::components(g)
        keep <- which(comp$membership == which.max(comp$csize))
        g <- igraph::induced_subgraph(g, keep)
    }
    adj <- methods::as(
        igraph::as_adjacency_matrix(g, sparse = TRUE), "CsparseMatrix")
    adj <- methods::as(adj, "generalMatrix")
    deg <- as.integer(Matrix::rowSums(adj))
    names(deg) <- igraph::V(g)$name
    methods::new("NFNetwork",
        graph = g, adjacency = methods::as(adj, "dMatrix"), degree = deg,
        bin = degreeBins(deg, minBinSize),
        minBinSize = as.integer(minBinSize), source = source)
}

#' Built-in STRING v12 human interactome
#'
#' @param evidence \code{"noKB"} (default) uses a combined score recomputed
#'   without STRING's curated-database and text-mining channels, so that the
#'   network does not encode the pathway annotations being tested;
#'   \code{"all"} uses the published combined score.
#' @param scoreThreshold minimum score (0--1000). Defaults: 400 for
#'   \code{"noKB"} and 700 for \code{"all"} (similar network sizes).
#' @param ... passed to \code{\link{buildNetwork}}.
#' @returns an \code{\linkS4class{NFNetwork}}.
#' @seealso \code{\link{STRING_v12_PPI}}
#' @export
#' @examples
#' stringNetwork()
stringNetwork <- function(evidence = c("noKB", "all"),
    scoreThreshold = NULL, ...) {
    evidence <- match.arg(evidence)
    env <- new.env(parent = emptyenv())
    utils::data("STRING_v12_PPI", package = "NetFunrankR", envir = env)
    ppi <- env$STRING_v12_PPI
    col <- if (evidence == "noKB") "nokb_score" else "combined_score"
    if (is.null(scoreThreshold))
        scoreThreshold <- if (evidence == "noKB") 400 else 700
    if (evidence == "noKB" && scoreThreshold < 400)
        warning("bundled noKB edges are only complete for scores >= 400")
    src <- sprintf("STRING v12 (%s evidence, score >= %d)",
        if (evidence == "noKB") "no database/text-mining" else "all",
        as.integer(scoreThreshold))
    buildNetwork(ppi[ppi[[col]] >= scoreThreshold, c("from", "to")],
        source = src, ...)
}

#' Degree binning (Guney et al. 2016)
#'
#' @param deg integer vector of node degrees.
#' @param minBinSize minimum number of nodes per bin.
#' @returns integer vector of bin labels (1 = lowest degree), in the order
#'   of \code{deg}.
#' @export
#' @examples
#' degreeBins(c(1, 1, 2, 2, 3, 10, 50), minBinSize = 2)
degreeBins <- function(deg, minBinSize = 100L) {
    tab <- table(deg)
    dvals <- as.integer(names(tab))
    counts <- as.integer(tab)
    binOfValue <- integer(length(dvals))
    b <- 1L
    acc <- 0L
    for (i in seq_along(dvals)) {
        binOfValue[i] <- b
        acc <- acc + counts[i]
        if (acc >= minBinSize && i < length(dvals)) {
            b <- b + 1L
            acc <- 0L
        }
    }
    last <- binOfValue == b
    if (b > 1L && sum(counts[last]) < minBinSize)
        binOfValue[last] <- b - 1L
    binOfValue[match(deg, dvals)]
}

#' @rdname NFNetwork-accessors
setMethod("nodeNames", "NFNetwork", function(x, ...) names(x@degree))

#' @rdname NFNetwork-accessors
setMethod("nodeDegree", "NFNetwork", function(x, ...) x@degree)

#' @rdname NFNetwork-accessors
setMethod("degreeBin", "NFNetwork", function(x, ...) x@bin)

#' @rdname NFNetwork-accessors
setMethod("adjacency", "NFNetwork", function(x, ...) x@adjacency)

#' @rdname NFNetwork-accessors
setMethod("nfGraph", "NFNetwork", function(x, ...) x@graph)

#' @rdname NFNetwork-accessors
#' @param object an \code{\linkS4class{NFNetwork}}.
setMethod("show", "NFNetwork", function(object) {
    cat("NFNetwork:", object@source, "\n")
    cat(sprintf("  %d nodes, %d edges (largest connected component)\n",
        length(object@degree), igraph::ecount(object@graph)))
    cat(sprintf("  %d degree bins (minimum %d nodes); median degree %g\n",
        length(unique(object@bin)), object@minBinSize,
        stats::median(object@degree)))
})

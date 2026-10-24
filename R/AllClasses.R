#' NFNetwork: an interactome prepared for network-proximity tests
#'
#' Holds the largest connected component of an undirected protein-protein
#' interaction network together with its sparse adjacency matrix, node
#' degrees and degree bins. Create objects with \code{\link{buildNetwork}}
#' or \code{\link{stringNetwork}}.
#'
#' @slot graph an \code{igraph} object (simple, undirected).
#' @slot adjacency sparse adjacency matrix (\code{dgCMatrix}).
#' @slot degree integer vector of node degrees, named by node.
#' @slot bin integer vector of degree-bin labels (1 = lowest degree).
#' @slot minBinSize minimum number of nodes per degree bin.
#' @slot source character description of the network source.
#'
#' @seealso \code{\link{buildNetwork}}, \code{\link{nodeNames}}
#' @exportClass NFNetwork
#' @importClassesFrom Matrix dgCMatrix
#' @examples
#' net <- stringNetwork()
#' net
#' head(nodeDegree(net))
setClass("NFNetwork",
    slots = c(
        graph = "ANY",
        adjacency = "dgCMatrix",
        degree = "integer",
        bin = "integer",
        minBinSize = "integer",
        source = "character"
    )
)

setValidity("NFNetwork", function(object) {
    n <- length(object@degree)
    msg <- character()
    if (!inherits(object@graph, "igraph"))
        msg <- c(msg, "'graph' must be an igraph object")
    if (nrow(object@adjacency) != n || ncol(object@adjacency) != n)
        msg <- c(msg, "'adjacency' must be a square matrix over all nodes")
    if (length(object@bin) != n)
        msg <- c(msg, "'bin' must have one entry per node")
    if (is.null(names(object@degree)) || anyNA(names(object@degree)))
        msg <- c(msg, "'degree' must be named by node")
    if (length(msg)) msg else TRUE
})

#' NFResult: results of a decomposed network-proximity enrichment analysis
#'
#' Returned by \code{\link{netEnrich}}. Use \code{\link{enrichResults}} to
#' obtain the per-gene-set table and \code{\link{evidenceClass}} to
#' attribute significant sets to overlap and/or network evidence.
#'
#' @slot results \code{\link[S4Vectors]{DataFrame}} with one row per tested
#'   gene set (see \code{\link{netEnrich}} for the columns).
#' @slot query character vector of query genes found in the network.
#' @slot unmatched character vector of query genes absent from the network.
#' @slot geneSets named list of the tested gene sets.
#' @slot params list of the analysis parameters.
#'
#' @exportClass NFResult
#' @importClassesFrom S4Vectors DataFrame
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' kegg <- loadPathways("KEGG")
#' res <- netEnrich(COPD_proteinCodingGenes, kegg, net)
#' res
setClass("NFResult",
    slots = c(
        results = "DataFrame",
        query = "character",
        unmatched = "character",
        geneSets = "list",
        params = "list"
    )
)

setValidity("NFResult", function(object) {
    need <- c("set", "p_overlap", "p_linker", "p_combined", "fdr_combined")
    miss <- setdiff(need, colnames(object@results))
    if (length(miss))
        return(paste("missing result columns:", paste(miss, collapse = ", ")))
    TRUE
})

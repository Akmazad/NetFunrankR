#' Accessors for NFNetwork objects
#'
#' @param x an \code{\linkS4class{NFNetwork}}.
#' @param ... not used.
#' @returns \code{nodeNames}: character vector of node (gene) names;
#'   \code{nodeDegree}: named integer vector of degrees;
#'   \code{degreeBin}: integer vector of degree-bin labels;
#'   \code{adjacency}: sparse adjacency matrix;
#'   \code{nfGraph}: the underlying \code{igraph} object.
#' @name NFNetwork-accessors
#' @examples
#' net <- stringNetwork()
#' length(nodeNames(net))
#' table(degreeBin(net))[1:5]
NULL

#' @rdname NFNetwork-accessors
#' @export
setGeneric("nodeNames", function(x, ...) standardGeneric("nodeNames"))

#' @rdname NFNetwork-accessors
#' @export
setGeneric("nodeDegree", function(x, ...) standardGeneric("nodeDegree"))

#' @rdname NFNetwork-accessors
#' @export
setGeneric("degreeBin", function(x, ...) standardGeneric("degreeBin"))

#' @rdname NFNetwork-accessors
#' @export
setGeneric("adjacency", function(x, ...) standardGeneric("adjacency"))

#' @rdname NFNetwork-accessors
#' @export
setGeneric("nfGraph", function(x, ...) standardGeneric("nfGraph"))

#' Accessors for NFResult objects
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param ... not used.
#' @returns \code{enrichResults}: a \code{DataFrame} with one row per gene
#'   set; \code{queryGenes}: query genes found in the network;
#'   \code{unmatchedGenes}: query genes not in the network;
#'   \code{nfParams}: list of analysis parameters;
#'   \code{testedGeneSets}: named list of the tested gene sets.
#' @name NFResult-accessors
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' head(enrichResults(res), 3)
#' length(queryGenes(res))
NULL

#' @rdname NFResult-accessors
#' @export
setGeneric("enrichResults", function(x, ...) standardGeneric("enrichResults"))

#' @rdname NFResult-accessors
#' @export
setGeneric("queryGenes", function(x, ...) standardGeneric("queryGenes"))

#' @rdname NFResult-accessors
#' @export
setGeneric("unmatchedGenes",
    function(x, ...) standardGeneric("unmatchedGenes"))

#' @rdname NFResult-accessors
#' @export
setGeneric("nfParams", function(x, ...) standardGeneric("nfParams"))

#' @rdname NFResult-accessors
#' @export
setGeneric("testedGeneSets",
    function(x, ...) standardGeneric("testedGeneSets"))

#' Evidence class of each gene set
#'
#' Gene sets significant by the combined test (FDR < \code{alpha}) are
#' labelled by which component reaches nominal significance
#' (p < \code{alpha}): the degree-aware overlap test, the beyond-overlap
#' linker test, both, or neither alone ("joint").
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param alpha significance threshold.
#' @param ... not used.
#' @returns factor with levels "overlap + network", "overlap only",
#'   "network only", "joint" and "none", one entry per gene set.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' table(evidenceClass(res))
setGeneric("evidenceClass",
    function(x, alpha = 0.05, ...) standardGeneric("evidenceClass"))

#' NetFunrankR: decomposed network-proximity functional enrichment
#'
#' Tests whether a gene list is closer to each pathway in a protein
#' interaction network than degree-matched random genes and decomposes the
#' evidence into a degree-aware overlap test and beyond-overlap (linker and
#' distance) tests of the non-member genes. All null distributions are
#' exact; no permutations are needed.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{stringNetwork}}, \code{\link{buildNetwork}}:
#'     prepare an interactome (\code{\linkS4class{NFNetwork}}).
#'   \item \code{\link{loadPathways}}, \code{\link{asGeneSetList}},
#'     \code{\link{readGMT}}: gene-set libraries.
#'   \item \code{\link{netEnrich}}, \code{\link{netEnrichBatch}}: the test
#'     (\code{\linkS4class{NFResult}}).
#'   \item \code{\link{evidenceClass}}, \code{\link{plotEvidenceMap}},
#'     \code{\link{plotTopSets}}, \code{\link{plotLinkerNetwork}}:
#'     interpretation.
#'   \item \code{\link{runNetFunrankRApp}}: interactive dashboard.
#' }
#'
#' @keywords internal
#' @useDynLib NetFunrankR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @importFrom methods new is as setClass setGeneric setMethod setValidity
#' @importFrom S4Vectors DataFrame
"_PACKAGE"

## variables evaluated inside mirai expressions in the dashboard
utils::globalVariables(c("lp", "compute", "setConfig", "cfg"))

## bundled example genes, for the dashboard
.exampleGenes <- function() {
    env <- new.env(parent = emptyenv())
    utils::data("COPD_proteinCodingGenes", package = "NetFunrankR",
        envir = env)
    env$COPD_proteinCodingGenes
}

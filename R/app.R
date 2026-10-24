.appPackages <- c("shiny", "bslib", "bsicons", "plotly", "reactable",
    "visNetwork", "shinyvalidate", "htmltools")

#' Interactive NetFunrankR dashboard
#'
#' Launches a Shiny application for interactive analysis: paste or upload a
#' gene list, choose a gene-set library and interactome, run the decomposed
#' network-proximity test, explore the evidence map, results table and
#' linker networks, and download the results, an HTML report and the R code
#' that reproduces the analysis.
#'
#' The dashboard is not tied to the bundled human data: supply your own
#' gene-set libraries, interactomes (any organism or identifier type) and
#' example list to serve a customised instance, e.g. for a lab or a
#' consortium. Computation runs asynchronously in background worker
#' processes when the \pkg{mirai} package is installed, so the interface
#' stays responsive and several users can be served at once.
#'
#' @param geneSets named list of gene-set libraries offered in the app. Each
#'   element is anything accepted by \code{\link{asGeneSetList}} (named list,
#'   two-column table, GMT file path, \code{GeneSetCollection}) or a function
#'   without arguments returning one (loaded on first use). Default: the
#'   bundled KEGG, Reactome and WikiPathways libraries. Users can always
#'   upload their own GMT file in addition.
#' @param networks named list of interactomes offered in the app. Each
#'   element is an \code{\linkS4class{NFNetwork}}, an edge list accepted by
#'   \code{\link{buildNetwork}}, or a function returning one. Default: the
#'   bundled STRING v12 networks. Users can always upload an edge list.
#' @param exampleGenes optional character vector offered as a one-click
#'   example (default: the bundled COPD genes; \code{character(0)} hides the
#'   example).
#' @param exampleLabel short name of the example shown in the interface.
#' @param geneLabel label of the gene input box (identifier type expected).
#' @param workers number of background worker processes (with \pkg{mirai}).
#' @param maxUploadMB maximum upload size in megabytes.
#' @param ... passed to \code{\link[shiny]{runApp}} (e.g. \code{host},
#'   \code{port}, \code{launch.browser}).
#' @returns \code{runNetFunrankRApp}: called for its side effect (runs the
#'   app); \code{netFunrankRApp}: a \code{shiny.appobj}, e.g. for deployment
#'   on Shiny Server, Posit Connect, shinyapps.io or in a Docker container.
#' @seealso \code{vignette("NetFunrankR")} and the \code{deploy/} directory
#'   of the source repository for deployment recipes.
#' @export
#' @examples
#' if (interactive()) {
#'     runNetFunrankRApp()
#'
#'     ## a customised instance: own library and network, no example
#'     myNet <- buildNetwork(read.csv("my_edges.csv"))
#'     runNetFunrankRApp(geneSets = list("My pathways" = "my_sets.gmt"),
#'         networks = list("My interactome" = myNet),
#'         exampleGenes = character(0), geneLabel = "Ensembl gene IDs")
#' }
runNetFunrankRApp <- function(geneSets = NULL, networks = NULL,
    exampleGenes = NULL, exampleLabel = "COPD",
    geneLabel = "Gene symbols (HGNC)", workers = 1L, maxUploadMB = 100,
    ...) {
    app <- netFunrankRApp(geneSets = geneSets, networks = networks,
        exampleGenes = exampleGenes, exampleLabel = exampleLabel,
        geneLabel = geneLabel, workers = workers, maxUploadMB = maxUploadMB)
    shiny::runApp(app, ...)
}

#' @rdname runNetFunrankRApp
#' @export
netFunrankRApp <- function(geneSets = NULL, networks = NULL,
    exampleGenes = NULL, exampleLabel = "COPD",
    geneLabel = "Gene symbols (HGNC)", workers = 1L, maxUploadMB = 100) {
    miss <- .appPackages[!vapply(.appPackages, requireNamespace,
        logical(1), quietly = TRUE)]
    if (length(miss))
        stop("The dashboard needs these packages: ",
            paste(miss, collapse = ", "), ".\nInstall them with ",
            "BiocManager::install(c(\"", paste(miss, collapse = "\", \""),
            "\"))", call. = FALSE)
    cfg <- .appConfig(geneSets, networks, exampleGenes, exampleLabel,
        geneLabel, workers, maxUploadMB)
    .setAppConfig(cfg)
    www <- system.file("app", "www", package = "NetFunrankR")
    if (nzchar(www))
        shiny::addResourcePath("nfr-assets", www)
    shiny::shinyApp(ui = .appUI(cfg),
        server = function(input, output, session)
            .appServer(input, output, session, cfg),
        onStart = function() .appOnStart(cfg))
}

#' Default libraries and networks of the dashboard
#'
#' The gene-set libraries and interactomes offered by
#' \code{\link{runNetFunrankRApp}} when none are supplied. Use them as a
#' starting point to add or remove entries for a customised dashboard.
#'
#' @returns list with elements \code{geneSets} and \code{networks}, each a
#'   named list of functions that load the resource on first use.
#' @export
#' @examples
#' d <- dashboardDefaults()
#' names(d$geneSets)
#' names(d$networks)
#' if (interactive()) {
#'     ## offer KEGG only, plus a custom GMT file
#'     runNetFunrankRApp(geneSets = c(d$geneSets["KEGG 2021"],
#'         list("My sets" = "my_sets.gmt")))
#' }
dashboardDefaults <- function() {
    list(geneSets = .defaultLibraries(), networks = .defaultNetworks())
}

## Start background workers and apply upload limits for this app.
.appOnStart <- function(cfg) {
    old <- options(shiny.maxRequestSize = cfg$maxUploadMB * 1024^2)
    shiny::onStop(function() options(old))
    if (.appAsyncAvailable() && !mirai::daemons_set()) {
        mirai::daemons(cfg$workers, dispatcher = cfg$workers > 1L)
        ## workers load the same installed package and app configuration
        mirai::everywhere({
            .libPaths(lp)
            setConfig(cfg)
        }, lp = .libPaths(), setConfig = .setAppConfig, cfg = cfg)
        shiny::onStop(function() mirai::daemons(0L))
    }
}

## mirai workers load the installed package. A development copy loaded with
## pkgload has no Meta/package.rds and is run synchronously instead.
.appAsyncAvailable <- function() {
    path <- getNamespaceInfo("NetFunrankR", "path")
    isTRUE(getOption("NetFunrankR.async", TRUE)) &&
        requireNamespace("mirai", quietly = TRUE) &&
        file.exists(file.path(path, "Meta", "package.rds"))
}

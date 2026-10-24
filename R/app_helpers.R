## Helper functions of the dashboard without Shiny dependencies (tested in
## tests/testthat/test-app.R).

.upload <- ".upload"   # choice value of "upload your own" entries

## Bundled defaults. Each entry is a function (loaded lazily) carrying the R
## code that reproduces it, used for the "Reproduce in R" panel.
.withCode <- function(f, code) structure(f, code = code)

.defaultLibraries <- function() {
    list(
        "KEGG 2021" = .withCode(function() loadPathways("KEGG"),
            "loadPathways(\"KEGG\")"),
        "Reactome 2024" = .withCode(function() loadPathways("Reactome"),
            "loadPathways(\"Reactome\")"),
        "WikiPathways 2024" = .withCode(
            function() loadPathways("WikiPathways"),
            "loadPathways(\"WikiPathways\")"))
}

.defaultNetworks <- function() {
    list(
        "STRING v12, no database/text mining, score >= 400 (recommended)" =
            .withCode(function() stringNetwork("noKB", 400),
                "stringNetwork(\"noKB\", 400)"),
        "STRING v12, no database/text mining, score >= 700" =
            .withCode(function() stringNetwork("noKB", 700),
                "stringNetwork(\"noKB\", 700)"),
        "STRING v12, all evidence, score >= 700" =
            .withCode(function() stringNetwork("all", 700),
                "stringNetwork(\"all\", 700)"))
}

## Validate and assemble the app configuration.
.appConfig <- function(geneSets, networks, exampleGenes, exampleLabel,
    geneLabel, workers, maxUploadMB) {
    defaults <- is.null(geneSets) && is.null(networks)
    if (is.null(geneSets)) geneSets <- .defaultLibraries()
    if (is.null(networks)) networks <- .defaultNetworks()
    for (x in list(geneSets, networks))
        if (!is.list(x) || !length(x) || is.null(names(x)) ||
            any(!nzchar(names(x))) || anyDuplicated(names(x)))
            stop("'geneSets' and 'networks' must be non-empty lists with ",
                "unique names")
    if (is.null(exampleGenes)) exampleGenes <- .exampleGenes()
    list(geneSets = geneSets, networks = networks,
        exampleGenes = unique(as.character(exampleGenes)),
        exampleLabel = exampleLabel, geneLabel = geneLabel,
        workers = max(1L, as.integer(workers)), maxUploadMB = maxUploadMB,
        defaults = defaults)
}

## The active configuration (set in the app process and in every worker).
.appState <- new.env(parent = emptyenv())
.setAppConfig <- function(cfg) {
    assign("cfg", cfg, envir = .appState)
    invisible(TRUE)
}
.getAppConfig <- function() {
    if (!exists("cfg", envir = .appState, inherits = FALSE))
        .setAppConfig(.appConfig(NULL, NULL, NULL, "COPD", "Gene symbols",
            1L, 100))
    get("cfg", envir = .appState, inherits = FALSE)
}

## Materialise a configured library or network (functions are called).
.resolve <- function(x) if (is.function(x)) x() else x
.resolveNetwork <- function(x) {
    x <- .resolve(x)
    if (methods::is(x, "NFNetwork")) x else buildNetwork(x)
}

## Split pasted or uploaded text into unique gene identifiers.
.parseGenes <- function(text) {
    if (is.null(text) || !length(text)) return(character())
    g <- unlist(strsplit(paste(text, collapse = "\n"), "[\\s,;|]+",
        perl = TRUE))
    g <- trimws(gsub("^[\"']|[\"']$", "", g))
    unique(g[nzchar(g)])
}

## Read genes from an uploaded txt/csv/tsv file (all tokens).
.readGeneFile <- function(path, name = path) {
    if (grepl("\\.(csv|tsv)$", name, ignore.case = TRUE)) {
        sep <- if (grepl("\\.tsv$", name, ignore.case = TRUE)) "\t" else ","
        df <- utils::read.table(path, sep = sep, header = FALSE,
            quote = "\"", fill = TRUE, stringsAsFactors = FALSE)
        return(.parseGenes(unlist(df, use.names = FALSE)))
    }
    .parseGenes(readLines(path, warn = FALSE))
}

## Cache of networks, libraries and gene-set distances (per process), so that
## repeated analyses only pay for the test itself.
.appCache <- new.env(parent = emptyenv())

.cached <- function(key, expr) {
    if (!exists(key, envir = .appCache, inherits = FALSE))
        assign(key, expr, envir = .appCache)
    get(key, envir = .appCache, inherits = FALSE)
}

.appNetworkFor <- function(choice, edges = NULL, hash = NULL) {
    if (identical(choice, .upload))
        return(.cached(paste0("net:upload:", hash),
            buildNetwork(edges, source = "uploaded edge list")))
    cfg <- .getAppConfig()
    if (is.null(cfg$networks[[choice]]))
        stop("unknown network: ", choice)
    .cached(paste0("net:", choice), .resolveNetwork(cfg$networks[[choice]]))
}

.appLibraryFor <- function(choice, uploaded = NULL, hash = NULL) {
    if (identical(choice, .upload))
        return(.cached(paste0("gs:upload:", hash), asGeneSetList(uploaded)))
    cfg <- .getAppConfig()
    if (is.null(cfg$geneSets[[choice]]))
        stop("unknown gene-set library: ", choice)
    .cached(paste0("gs:", choice),
        asGeneSetList(.resolve(cfg$geneSets[[choice]])))
}

## Run one analysis described by 'spec' (a plain list that can be sent to a
## background worker).
.appCompute <- function(spec) {
    t0 <- proc.time()[["elapsed"]]
    net <- .appNetworkFor(spec$network, spec$edges, spec$networkHash)
    gs <- .appLibraryFor(spec$library, spec$geneSets, spec$gmtHash)
    netKey <- if (identical(spec$network, .upload)) spec$networkHash else
        spec$network
    libKey <- if (identical(spec$library, .upload)) spec$gmtHash else
        spec$library
    D <- .cached(paste("D", netKey, libKey, sep = ":"),
        geneSetDistances(net, gs))
    res <- netEnrich(spec$genes, gs, net, distances = D,
        minSize = spec$minSize, maxSize = spec$maxSize,
        background = spec$background, beyond = spec$beyond)
    list(result = res, seconds = proc.time()[["elapsed"]] - t0,
        spec = spec[c("genes", "network", "library", "minSize", "maxSize",
            "alpha", "beyond")])
}

## R code reproducing an entry: its 'code' attribute or a placeholder.
.entryCode <- function(entry, upload, placeholder) {
    if (upload) return(placeholder)
    code <- attr(entry, "code")
    if (is.null(code)) placeholder else code
}

## R code that reproduces a dashboard analysis.
.rCode <- function(spec, genes, cfg = .getAppConfig()) {
    q <- strwrap(paste0("\"", genes, "\"", collapse = ", "), width = 72)
    netCall <- .entryCode(cfg$networks[[spec$network]],
        identical(spec$network, .upload),
        "buildNetwork(read.csv(\"my_edges.csv\"))  # your network")
    gsCall <- .entryCode(cfg$geneSets[[spec$library]],
        identical(spec$library, .upload),
        "readGMT(\"my_gene_sets.gmt\")  # your gene sets")
    c("library(NetFunrankR)", "",
        paste0("genes <- c(", paste(q, collapse = "\n    "), ")"),
        paste0("net <- ", netCall),
        paste0("geneSets <- ", gsCall),
        sprintf(paste0("res <- netEnrich(genes, geneSets, net, minSize = %d,",
            " maxSize = %s,\n    beyond = \"%s\")"), as.integer(spec$minSize),
            format(spec$maxSize), spec$beyond),
        "res", sprintf("table(evidenceClass(res, alpha = %s))", spec$alpha),
        sprintf("plotEvidenceMap(res, alpha = %s)", spec$alpha),
        "exportResults(res, \"netfunrankr_results.csv\")")
}

## Compact display of p-values
.fmtP <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 1e-3, formatC(p, format = "e",
        digits = 1), formatC(p, format = "f", digits = 3)))
}

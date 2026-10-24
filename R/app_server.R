.appServer <- function(input, output, session, cfg = .getAppConfig()) {
    st <- .serverInputs(input, output, session, cfg)
    st$task <- .serverTask(input, session, st, cfg)
    st <- .serverResults(input, output, session, st)
    .serverPlots(input, output, session, st)
    .serverTable(input, output, session, st)
    .serverExplorer(input, output, session, st)
    .serverDownloads(input, output, session, st)
}

## ---------------------------------------------------------------- inputs
.serverInputs <- function(input, output, session, cfg) {
    genes <- shiny::reactive({
        g <- .parseGenes(input$genes)
        f <- input$genesFile
        if (!is.null(f)) g <- unique(c(g, .readGeneFile(f$datapath, f$name)))
        g
    })
    example <- paste(cfg$exampleGenes, collapse = "\n")
    shiny::observeEvent(input$example,
        shiny::updateTextAreaInput(session, "genes", value = example))
    shiny::observeEvent(input$clear,
        shiny::updateTextAreaInput(session, "genes", value = ""))
    network <- shiny::reactive({
        shiny::req(input$network)
        if (identical(input$network, .upload)) {
            shiny::req(input$edgeFile)
            hash <- unname(tools::md5sum(input$edgeFile$datapath))
            return(.appNetworkFor(.upload,
                .readEdges(input$edgeFile$datapath), hash))
        }
        .appNetworkFor(input$network)
    })
    output$geneStatus <- shiny::renderUI(.geneStatusUI(genes(), network()))
    iv <- .validator(input, genes)
    list(genes = genes, network = network, iv = iv)
}

.readEdges <- function(path) {
    first <- readLines(path, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", first)) "\t" else if (grepl(",", first)) "," else ""
    df <- utils::read.table(path, sep = sep, header = FALSE, quote = "\"",
        stringsAsFactors = FALSE, fill = TRUE)
    header <- c("from", "to", "gene1", "gene2", "protein1", "protein2",
        "source", "target")
    if (nrow(df) && all(tolower(df[1, c(1L, 2L)]) %in% header))
        df <- df[-1, ]
    df[, c(1L, 2L)]
}

.geneStatusUI <- function(genes, network) {
    if (!length(genes)) return(NULL)
    found <- genes %in% nodeNames(network)
    miss <- genes[!found]
    htmltools::div(class = "d-flex flex-wrap gap-1 small",
        htmltools::span(class = "badge text-bg-secondary",
            sprintf("%d genes", length(genes))),
        htmltools::span(class = "badge text-bg-success",
            sprintf("%d in network", sum(found))),
        if (length(miss)) bslib::popover(
            htmltools::span(class = "badge text-bg-warning nfr-clickable",
                role = "button", tabindex = "0",
                sprintf("%d not found", length(miss))),
            title = "Not in the network",
            htmltools::p(class = "small mb-0", paste(utils::head(miss, 200),
                collapse = ", "))))
}

.validator <- function(input, genes) {
    iv <- shinyvalidate::InputValidator$new()
    iv$add_rule("genes", function(value) {
        if (length(genes()) < 2) "Enter or upload at least two genes."
    })
    iv$add_rule("alpha", shinyvalidate::sv_between(0.0001, 0.5))
    iv$add_rule("gmtFile", function(value) {
        if (identical(input$library, .upload) && is.null(value))
            "Upload a GMT file or choose a library."
    })
    iv$add_rule("edgeFile", function(value) {
        if (identical(input$network, .upload) && is.null(value))
            "Upload an edge list or choose a network."
    })
    iv
}

## ---------------------------------------------------------------- task
.serverTask <- function(input, session, st, cfg) {
    task <- shiny::ExtendedTask$new(function(spec) {
        if (.appAsyncAvailable())
            return(mirai::mirai(compute(spec), compute = .appCompute,
                spec = spec))
        promises::promise(function(resolve, reject) {
            tryCatch(resolve(.appCompute(spec)), error = reject)
        })
    })
    bslib::bind_task_button(task, "run")
    launch <- function(genes) task$invoke(.buildSpec(input, genes))
    shiny::observeEvent(input$run, {
        st$iv$enable()
        if (!st$iv$is_valid()) {
            .toast("Please check the highlighted inputs.", "warning")
            return()
        }
        launch(st$genes())
    })
    shiny::observeEvent(input$tryExample, {
        shiny::updateTextAreaInput(session, "genes",
            value = paste(cfg$exampleGenes, collapse = "\n"))
        launch(cfg$exampleGenes)
    })
    task
}

.buildSpec <- function(input, genes) {
    spec <- list(genes = genes, network = input$network,
        library = input$library, minSize = input$sizeRange[1],
        maxSize = input$sizeRange[2], alpha = input$alpha,
        beyond = input$beyond, background = NULL)
    if (identical(input$library, .upload)) {
        spec$geneSets <- readGMT(input$gmtFile$datapath)
        spec$gmtHash <- unname(tools::md5sum(input$gmtFile$datapath))
    }
    if (identical(input$network, .upload)) {
        spec$edges <- .readEdges(input$edgeFile$datapath)
        spec$networkHash <- unname(tools::md5sum(input$edgeFile$datapath))
    }
    if (!is.null(input$bgFile))
        spec$background <- .readGeneFile(input$bgFile$datapath,
            input$bgFile$name)
    spec
}

.toast <- function(text, type = "primary") {
    bslib::show_toast(bslib::toast(text, type = type, duration_s = 6,
        position = "bottom-right"))
}

## ---------------------------------------------------------------- results
.serverResults <- function(input, output, session, st) {
    last <- shiny::reactiveVal(NULL)
    shiny::observe({
        status <- st$task$status()
        if (status == "success") {
            out <- st$task$result()
            last(out)
            shiny::updateSelectizeInput(session, "pathway", server = TRUE,
                choices = enrichResults(out$result)$set,
                selected = .defaultPathway(out$result, input$alpha))
            .toast(sprintf("Analysis finished in %.1f s.", out$seconds),
                "success")
        } else if (status == "error") {
            msg <- tryCatch({
                st$task$result()
                ""
            }, error = conditionMessage)
            .toast(paste("Analysis failed:", msg), "danger")
        }
    })
    st$out <- last
    st$res <- shiny::reactive({
        shiny::req(last())
        last()$result
    })
    st$alpha <- shiny::reactive({
        a <- input$alpha
        if (is.numeric(a) && !is.na(a) && a > 0 && a < 1) a else 0.05
    })
    output$hasResults <- shiny::reactive(!is.null(last()))
    shiny::outputOptions(output, "hasResults", suspendWhenHidden = FALSE)
    .valueBoxes(output, st)
    output$runInfo <- shiny::renderUI({
        o <- shiny::req(last())
        htmltools::p(class = "small text-body-secondary mt-2 mb-0",
            bsicons::bs_icon("clock-history"), sprintf("Last run: %.1f s",
                o$seconds), "\u00b7", o$spec$library, "\u00b7",
            nfParams(o$result)$network)
    })
    st
}

.defaultPathway <- function(res, alpha) {
    r <- enrichResults(res)
    ev <- evidenceClass(res, alpha)
    pick <- which(ev %in% c("network only", "overlap + network"))
    r$set[if (length(pick)) pick[1] else 1L]
}

## Reactives are read into local variables before S4 dispatch: methods
## turns shiny's silent req() condition into a visible error otherwise.
.valueBoxes <- function(output, st) {
    output$vbGenes <- shiny::renderText({
        r <- st$res()
        length(queryGenes(r))
    })
    output$vbGenesSub <- shiny::renderUI({
        r <- st$res()
        u <- length(unmatchedGenes(r))
        htmltools::p(if (u) sprintf("%d not in network", u) else
            "all genes matched")
    })
    output$vbSets <- shiny::renderText({
        r <- st$res()
        nrow(enrichResults(r))
    })
    output$vbSetsSub <- shiny::renderUI({
        r <- st$res()
        p <- nfParams(r)
        htmltools::p(sprintf("sizes %d\u2013%s", as.integer(p$minSize),
            format(p$maxSize)))
    })
    output$vbSig <- shiny::renderText({
        r <- st$res()
        sum(enrichResults(r)$fdr_combined < st$alpha())
    })
    output$vbSigSub <- shiny::renderUI({
        r <- st$res()
        htmltools::p(sprintf("FDR < %s; ORA finds %d", st$alpha(),
            sum(enrichResults(r)$fdr_hyper < st$alpha())))
    })
    output$vbNet <- shiny::renderText({
        r <- st$res()
        sum(evidenceClass(r, st$alpha()) == "network only")
    })
    output$vbNetSub <- shiny::renderUI(htmltools::p(
        "found through linker genes"))
}

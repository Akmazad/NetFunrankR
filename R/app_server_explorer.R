## ---------------------------------------------------------------- explorer
.serverExplorer <- function(input, output, session, st) {
    row <- shiny::reactive({
        res <- st$res()
        r <- enrichResults(res)
        shiny::req(input$pathway %in% r$set)
        i <- match(input$pathway, r$set)
        list(r = r[i, ], ev = as.character(evidenceClass(res,
            st$alpha())[i]))
    })
    output$pathwayStats <- shiny::renderUI(.pathwayStatsUI(row()))
    output$linkerNet <- visNetwork::renderVisNetwork({
        res <- st$res()
        shiny::req(input$pathway %in% names(testedGeneSets(res)))
        .linkerVis(res, input$pathway, st$network(), input$dark == "dark")
    })
    output$nullPlot <- plotly::renderPlotly({
        res <- st$res()
        shiny::req(input$pathway %in% names(testedGeneSets(res)))
        nd <- nullDistribution(queryGenes(res),
            testedGeneSets(res)[input$pathway], st$network(), input$pathway)
        .nullPlot(nd, input$dark == "dark")
    })
}

.pathwayStatsUI <- function(x) {
    r <- x$r
    chips <- function(g) {
        g <- strsplit(g, "/", fixed = TRUE)[[1]]
        if (!length(g)) return(htmltools::span(class = "text-body-secondary",
            "none"))
        lapply(g, function(z) htmltools::span(class = "nfr-chip", z))
    }
    stat <- function(label, value) htmltools::div(class = "nfr-stat",
        htmltools::span(class = "text-body-secondary", label),
        htmltools::span(class = "fw-semibold", value))
    htmltools::tagList(
        htmltools::span(class = "nfr-pill mb-2",
            style = sprintf("--pill: %s;", .evidenceCols[[x$ev]]), x$ev),
        stat("Genes in network", r$set_in_net),
        stat("Overlap (expected)", sprintf("%d (%.1f)", r$overlap,
            r$overlap_expected)),
        stat("Linker genes (expected)", sprintf("%d (%.1f)", r$n_linker,
            r$linker_expected)),
        stat("p overlap", .fmtP(r$p_overlap)),
        stat("p linker", .fmtP(r$p_linker)),
        stat("FDR combined", .fmtP(r$fdr_combined)),
        stat("FDR ORA", .fmtP(r$fdr_hyper)),
        htmltools::h6(class = "mt-3", "Query genes in the set"),
        htmltools::div(class = "nfr-chips", chips(r$overlap_genes)),
        htmltools::h6(class = "mt-3", "Linker genes"),
        htmltools::div(class = "nfr-chips", chips(r$linker_genes)))
}

.linkerVis <- function(res, set, network, dark) {
    shiny::req(set)
    sub <- linkerSubnetwork(res, set, network, maxMembers = 40L)
    cols <- .modeCols(dark)
    roleCol <- c("Query gene in set" = cols[["overlap + network"]],
        "Query linker gene" = cols[["network only"]],
        "Set member" = cols[["member"]])
    vn <- igraph::V(sub$graph)$name
    role <- sub$role[vn]
    nodes <- data.frame(id = vn, label = ifelse(role == "Set member", "", vn),
        title = paste0("<b>", vn, "</b><br>", role), group = role,
        color = unname(roleCol[role]), size = ifelse(role == "Set member",
            14, 24), font.color = cols[["ink"]], font.size = 18)
    el <- igraph::as_edgelist(sub$graph)
    edges <- data.frame(from = el[, 1], to = el[, 2], color = cols[["grid"]],
        width = 1.2)
    ## a stable layout computed in R; no physics simulation in the browser
    lay <- .stableLayout(sub$graph)
    lay <- apply(lay, 2, function(v) (v - min(v)) / max(diff(range(v)), 1e-9))
    nodes$x <- (lay[, 1] - 0.5) * 700
    nodes$y <- (lay[, 2] - 0.5) * 500
    v <- visNetwork::visNetwork(nodes, edges, width = "100%")
    v <- visNetwork::visNodes(v, shape = "dot", borderWidth = 1,
        physics = FALSE)
    v <- visNetwork::visEdges(v, physics = FALSE, smooth = FALSE)
    v <- visNetwork::visOptions(v, highlightNearest = list(enabled = TRUE,
        degree = 1, hover = TRUE))
    visNetwork::visInteraction(v, tooltipDelay = 50, zoomView = TRUE)
}

.nullPlot <- function(nd, dark) {
    p <- nd$S
    df <- data.frame(d = (seq_along(p) - 1) / nd$n, prob = p)
    df <- df[df$prob > 1e-6 * max(df$prob), ]
    obs <- nd$S_obs / nd$n
    cols <- .modeCols(dark)
    g <- plotly::plot_ly(df, x = ~d, y = ~prob, type = "bar",
        marker = list(color = cols[["member"]]), hoverinfo = "x+y",
        name = "null")
    g <- plotly::layout(g, bargap = 0.05, showlegend = FALSE,
        shapes = list(list(type = "line", x0 = obs, x1 = obs, y0 = 0,
            y1 = 1, yref = "paper", line = list(color =
                cols[["network only"]], width = 3))),
        annotations = list(list(x = obs, y = 1, yref = "paper",
            text = sprintf("observed %.2f", obs), showarrow = FALSE,
            xanchor = "left", font = list(color = cols[["network only"]]))))
    .plotlyStyle(g, dark, "Mean distance of query genes to the set (null)",
        "Probability")
}

## ---------------------------------------------------------------- downloads
.serverDownloads <- function(input, output, session, st) {
    code <- shiny::reactive({
        o <- shiny::req(st$out())
        .rCode(o$spec, o$spec$genes)
    })
    output$rcode <- shiny::renderText(paste(code(), collapse = "\n"))
    output$dlCsv <- shiny::downloadHandler(
        filename = function() sprintf("netfunrankr_%s.csv", Sys.Date()),
        content = function(file) {
            res <- st$res()
            exportResults(res, file, st$alpha())
        })
    output$dlCode <- shiny::downloadHandler(
        filename = function() "netfunrankr_analysis.R",
        content = function(file) writeLines(code(), file))
    output$dlReport <- shiny::downloadHandler(
        filename = function() sprintf("netfunrankr_report_%s.html",
            Sys.Date()),
        content = function(file) {
            shiny::withProgress(message = "Rendering report\u2026", {
                res <- st$res()
                .renderReport(res, st$network(), st$alpha(), code(), file)
            })
        })
}

.renderReport <- function(res, network, alpha, code, file) {
    if (!requireNamespace("rmarkdown", quietly = TRUE))
        stop("install 'rmarkdown' to create reports")
    tmp <- tempfile("nfr-report")
    dir.create(tmp)
    rmd <- file.path(tmp, "report.Rmd")
    file.copy(system.file("app", "report.Rmd", package = "NetFunrankR"), rmd)
    out <- rmarkdown::render(rmd, output_dir = tmp, quiet = TRUE,
        params = list(result = res, network = network, alpha = alpha,
            code = code), envir = new.env(parent = globalenv()))
    file.copy(out, file, overwrite = TRUE)
}

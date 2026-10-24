## Colours for light and dark mode (the dark steps of the same hues)
.modeCols <- function(dark) {
    if (isTRUE(dark))
        c("overlap + network" = "#199e70", "overlap only" = "#3987e5",
            "network only" = "#d95926", "joint" = "#9085e9",
            "none" = "#5d5c58", ink = "#f0efec", grid = "#4a4945",
            member = "#7e9aab")
    else c(.evidenceCols, ink = "#0b0b0b", grid = "#e8e7e3",
        member = "#9fb6c3")
}

.plotlyStyle <- function(p, dark, xlab, ylab, source = NULL) {
    cols <- .modeCols(dark)
    ax <- function(t) list(title = t, gridcolor = cols[["grid"]],
        zerolinecolor = cols[["grid"]], color = cols[["ink"]])
    p <- plotly::layout(p, xaxis = ax(xlab), yaxis = ax(ylab),
        paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
        font = list(color = cols[["ink"]],
            family = "Inter, system-ui, sans-serif"),
        legend = list(orientation = "h", y = -0.18),
        hoverlabel = list(align = "left"), margin = list(t = 20))
    plotly::config(p, displaylogo = FALSE,
        toImageButtonOptions = list(format = "svg", filename = "netfunrankr"))
}

.hoverText <- function(df) {
    sprintf(paste0("<b>%s</b><br>Evidence: %s<br>Overlap: %d (expected %.1f)",
        "<br>Linker genes: %d (expected %.1f)<br>FDR combined: %s",
        "<br>FDR ORA: %s"), df$set, df$evidence, df$overlap,
        df$overlap_expected, df$n_linker, df$linker_expected,
        .fmtP(df$fdr_combined), .fmtP(df$fdr_hyper))
}

## ---------------------------------------------------------------- plots
.serverPlots <- function(input, output, session, st) {
    output$evidenceMap <- plotly::renderPlotly({
        res <- st$res()
        df <- as.data.frame(res, alpha = st$alpha())
        yCol <- if (nfParams(res)$beyond == "distance") "p_beyond" else
            "p_linker"
        df$x <- -log10(df$p_overlap)
        df$y <- -log10(df[[yCol]])
        df$evidence <- factor(df$evidence, .evidenceLevels)
        df$sz <- 7 + 16 * sqrt(df$set_in_net / max(df$set_in_net))
        cols <- .modeCols(input$dark == "dark")
        top <- utils::head(df[df$evidence != "none", ], 5)
        p <- plotly::plot_ly(df, x = ~x, y = ~y, color = ~evidence,
            colors = cols[.evidenceLevels], type = "scatter",
            mode = "markers", text = .hoverText(df), hoverinfo = "text",
            customdata = ~set, source = "evmap",
            marker = list(size = ~sz, opacity = 0.85,
                line = list(width = 1, color = "rgba(255,255,255,0.7)")))
        if (nrow(top))
            p <- plotly::add_annotations(p, data = top, x = ~x, y = ~y,
                text = ~set, showarrow = TRUE, arrowhead = 0, ax = 28,
                ay = rep(c(-34, 30, -58, 52, -82), length.out = nrow(top)),
                arrowcolor = cols[["grid"]], bgcolor = "rgba(0,0,0,0)",
                font = list(size = 11, color = cols[["ink"]]), inherit = FALSE)
        p <- plotly::event_register(p, "plotly_click")
        .plotlyStyle(p, input$dark == "dark",
            "\u2212log\u2081\u2080 p overlap (degree-aware)",
            paste("\u2212log\u2081\u2080 p", sub("p_", "", yCol),
                "(beyond overlap)"))
    })
    ## observe plotly's raw input so no event is read before the plot exists
    shiny::observeEvent(input[["plotly_click-evmap"]], {
        cd <- plotly::event_data("plotly_click", source = "evmap")$customdata
        if (length(cd) && !is.na(cd[[1]])) .selectPathway(session, cd[[1]])
    })
    output$topSets <- plotly::renderPlotly({
        res <- st$res()
        .topSetsPlot(res, st$alpha(), input$dark == "dark")
    })
    output$oraCompare <- plotly::renderPlotly({
        res <- st$res()
        .oraPlot(res, st$alpha(), input$dark == "dark")
    })
}

.selectPathway <- function(session, set) {
    shiny::updateSelectizeInput(session, "pathway", selected = set)
    bslib::nav_select("tabs", "explorer", session = session)
}

.topSetsPlot <- function(res, alpha, dark) {
    df <- utils::head(as.data.frame(res, alpha = alpha), 25)
    df$set <- factor(df$set, rev(df$set))
    df$lp <- -log10(df$p_combined)
    cols <- .modeCols(dark)
    p <- plotly::plot_ly(df, y = ~set, hoverinfo = "text",
        text = .hoverText(df))
    p <- plotly::add_segments(p, x = 0, xend = ~lp, yend = ~set,
        line = list(color = cols[["grid"]], width = 2), showlegend = FALSE)
    p <- plotly::add_markers(p, x = ~lp, color = ~factor(evidence,
        .evidenceLevels), colors = cols[.evidenceLevels],
        marker = list(size = 12))
    .plotlyStyle(p, dark, "\u2212log\u2081\u2080 p combined", "")
}

.oraPlot <- function(res, alpha, dark) {
    df <- as.data.frame(res, alpha = alpha)
    ora <- df$fdr_hyper < alpha
    nfr <- df$fdr_combined < alpha
    df$group <- factor(ifelse(ora & nfr, "Both", ifelse(nfr,
        "NetFunrankR only", ifelse(ora, "ORA only", "Neither"))),
        c("Both", "NetFunrankR only", "ORA only", "Neither"))
    cols <- .modeCols(dark)
    pal <- c(Both = cols[["overlap + network"]],
        `NetFunrankR only` = cols[["network only"]],
        `ORA only` = cols[["overlap only"]], Neither = cols[["none"]])
    lim <- max(-log10(c(df$p_hyper, df$p_combined)), 1)
    p <- plotly::plot_ly(df, x = ~ -log10(p_hyper), y = ~ -log10(p_combined),
        color = ~group, colors = pal, type = "scatter", mode = "markers",
        text = .hoverText(df), hoverinfo = "text", marker = list(size = 9,
            opacity = 0.8))
    p <- plotly::layout(p, shapes = list(list(type = "line", x0 = 0, y0 = 0,
        x1 = lim, y1 = lim, line = list(dash = "dot",
            color = cols[["grid"]]))))
    .plotlyStyle(p, dark, "\u2212log\u2081\u2080 p hypergeometric ORA",
        "\u2212log\u2081\u2080 p NetFunrankR (combined)")
}

## ---------------------------------------------------------------- table
.serverTable <- function(input, output, session, st) {
    output$table <- reactable::renderReactable({
        res <- st$res()
        .resultsTable(as.data.frame(res, alpha = st$alpha()))
    })
    shiny::observeEvent(reactable::getReactableState("table", "selected"), {
        i <- reactable::getReactableState("table", "selected")
        res <- st$res()
        if (length(i)) .selectPathway(session, enrichResults(res)$set[i])
    })
}

.resultsTable <- function(df) {
    keep <- c("set", "evidence", "overlap", "n_linker", "set_in_net",
        "p_overlap", "p_linker", "fdr_combined", "fdr_hyper")
    d <- df[, keep]
    pcol <- function(name) reactable::colDef(name = name, minWidth = 90,
        cell = function(v) .fmtP(v), align = "right")
    reactable::reactable(d, searchable = TRUE, filterable = FALSE,
        highlight = TRUE, compact = TRUE, defaultPageSize = 15,
        selection = "single", onClick = "select", striped = FALSE,
        theme = .reactableTheme(),
        columns = list(
            set = reactable::colDef(name = "Gene set", minWidth = 260),
            evidence = reactable::colDef(name = "Evidence", minWidth = 150,
                cell = function(v) htmltools::span(class = "nfr-pill",
                    style = sprintf("--pill: %s;", .evidenceCols[[v]]), v)),
            overlap = reactable::colDef(name = "Overlap", align = "right",
                cell = function(v, i) sprintf("%d (%.1f)", v,
                    df$overlap_expected[i])),
            n_linker = reactable::colDef(name = "Linkers", align = "right",
                cell = function(v, i) sprintf("%d (%.1f)", v,
                    df$linker_expected[i])),
            set_in_net = reactable::colDef(name = "Size", align = "right",
                maxWidth = 70),
            p_overlap = pcol("p overlap"), p_linker = pcol("p linker"),
            fdr_combined = pcol("FDR combined"),
            fdr_hyper = pcol("FDR ORA")))
}

.reactableTheme <- function() {
    reactable::reactableTheme(color = "var(--bs-body-color)",
        backgroundColor = "transparent",
        borderColor = "var(--bs-border-color)",
        highlightColor = "var(--bs-tertiary-bg)",
        inputStyle = list(backgroundColor = "var(--bs-body-bg)",
            color = "var(--bs-body-color)"),
        selectStyle = list(backgroundColor = "var(--bs-body-bg)"),
        pageButtonStyle = list(color = "var(--bs-body-color)"))
}

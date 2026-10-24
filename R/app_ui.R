.appTheme <- function() {
    bslib::bs_theme(version = 5,
        primary = "#2a78d6", success = "#1baf7a", warning = "#eda100",
        danger = "#e34948", info = "#6250d6",
        base_font = bslib::font_collection("Inter", "system-ui",
            "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue",
            "Arial", "sans-serif"),
        code_font = bslib::font_collection("JetBrains Mono", "Consolas",
            "SFMono-Regular", "Menlo", "monospace"),
        "border-radius" = "0.6rem", "navbar-light-bg" = "#ffffff")
}

.appUI <- function(cfg = .getAppConfig()) {
    bslib::page_navbar(
        id = "nav", title = .brand(), theme = .appTheme(),
        window_title = "NetFunrankR", fillable = "analyze",
        header = htmltools::tags$head(
            htmltools::tags$link(rel = "stylesheet",
                href = "nfr-assets/styles.css"),
            htmltools::tags$script(src = "nfr-assets/app.js")),
        bslib::nav_panel("Analyze", value = "analyze",
            icon = bsicons::bs_icon("bullseye"),
            bslib::layout_sidebar(sidebar = .sidebarUI(cfg), fillable = TRUE,
                .welcomeUI(cfg), .resultsUI())),
        bslib::nav_panel("How it works", value = "guide",
            icon = bsicons::bs_icon("lightbulb"), .guideUI()),
        bslib::nav_panel("About", value = "about",
            icon = bsicons::bs_icon("info-circle"), .aboutUI(cfg)),
        bslib::nav_spacer(),
        bslib::nav_item(htmltools::tags$a(class = "nav-link",
            href = .projectURL(), target = "_blank", rel = "noopener",
            `aria-label` = "NetFunrankR source code",
            bsicons::bs_icon("github"))),
        bslib::nav_item(bslib::input_dark_mode(id = "dark", mode = "light"))
    )
}

## Project URL taken from the package DESCRIPTION.
.projectURL <- function() {
    url <- utils::packageDescription("NetFunrankR")$URL
    if (is.null(url)) "https://bioconductor.org" else
        trimws(strsplit(url, ",")[[1]][1])
}

.brand <- function() {
    htmltools::tags$span(class = "nfr-brand",
        htmltools::tags$img(src = "nfr-assets/logo.svg", height = 28,
            alt = ""),
        htmltools::tags$span("NetFunrankR"))
}

.sidebarUI <- function(cfg) {
    bslib::sidebar(width = 360, open = list(desktop = "open",
        mobile = "closed"),
        bslib::accordion(id = "setup", multiple = TRUE,
            open = c("genes", "library", "network"),
            .genesPanel(cfg), .libraryPanel(cfg), .networkPanel(cfg),
            .advancedPanel()),
        bslib::input_task_button("run", "Run analysis",
            icon = bsicons::bs_icon("play-fill"),
            label_busy = "Analysing...", class = "btn-primary w-100"),
        shiny::uiOutput("runInfo"))
}

.hasExample <- function(cfg) length(cfg$exampleGenes) >= 2L

.genesPanel <- function(cfg) {
    ex <- utils::head(cfg$exampleGenes, 4)
    bslib::accordion_panel("Query genes", value = "genes",
        icon = bsicons::bs_icon("list-ul"),
        shiny::textAreaInput("genes", cfg$geneLabel, rows = 7,
            width = "100%", resize = "vertical",
            placeholder = paste0("One per line or separated by spaces or ",
                "commas", if (length(ex)) paste0(", e.g.\n",
                    paste(ex, collapse = " ")))),
        shiny::fileInput("genesFile", NULL, buttonLabel = "Upload...",
            placeholder = "txt, csv or tsv", width = "100%",
            accept = c(".txt", ".csv", ".tsv")),
        htmltools::div(class = "d-flex justify-content-between mb-2",
            if (.hasExample(cfg)) shiny::actionLink("example",
                htmltools::tagList(bsicons::bs_icon("stars"),
                    paste("Load", cfg$exampleLabel, "example"))),
            shiny::actionLink("clear", htmltools::tagList(
                bsicons::bs_icon("x-circle"), "Clear"))),
        shiny::uiOutput("geneStatus"))
}

.choices <- function(entries, uploadLabel) {
    stats::setNames(c(names(entries), .upload),
        c(names(entries), uploadLabel))
}

.libraryPanel <- function(cfg) {
    bslib::accordion_panel("Gene-set library", value = "library",
        icon = bsicons::bs_icon("collection"),
        shiny::selectInput("library", NULL,
            .choices(cfg$geneSets, "Upload a GMT file"), width = "100%"),
        shiny::conditionalPanel(sprintf("input.library == '%s'", .upload),
            shiny::fileInput("gmtFile", "GMT file", accept = ".gmt",
                width = "100%")))
}

.networkPanel <- function(cfg) {
    bslib::accordion_panel("Interactome", value = "network",
        icon = bsicons::bs_icon("diagram-3"),
        shiny::selectInput("network", NULL,
            .choices(cfg$networks, "Upload an edge list"), width = "100%"),
        shiny::conditionalPanel(sprintf("input.network == '%s'", .upload),
            shiny::fileInput("edgeFile", "Edge list (first two columns)",
                accept = c(".csv", ".tsv", ".txt"), width = "100%")),
        if (cfg$defaults) htmltools::p(class = "small text-body-secondary mb-0",
            "The recommended network excludes STRING's curated-database and",
            "text-mining evidence, so it does not encode the pathway",
            "annotations being tested."))
}

.advancedPanel <- function() {
    bslib::accordion_panel("Advanced settings", value = "advanced",
        icon = bsicons::bs_icon("sliders"),
        shiny::sliderInput("sizeRange", "Gene-set size (genes in network)",
            min = 1, max = 1000, value = c(5, 500), step = 1, width = "100%"),
        shiny::numericInput("alpha", "FDR threshold", 0.05, min = 0.001,
            max = 0.25, step = 0.01, width = "100%"),
        shiny::radioButtons("beyond", .labelWithHelp("Beyond-overlap test",
            paste("Linker: number of query genes that interact directly",
                "with the gene set without being annotated to it.",
                "Distance: their mean shortest-path distance.")),
            c("Linker genes (recommended)" = "linker",
                "Mean distance" = "distance")),
        shiny::fileInput("bgFile", .labelWithHelp("Background genes",
            paste("Optional: genes that could have been selected, e.g. all",
                "genes measured in an expression experiment.")),
            accept = c(".txt", ".csv", ".tsv"), width = "100%"))
}

.labelWithHelp <- function(label, help) {
    htmltools::tagList(label, bslib::tooltip(
        bsicons::bs_icon("question-circle", title = help,
            class = "ms-1 text-body-secondary"), help))
}

.welcomeUI <- function(cfg) {
    libs <- paste(utils::head(names(cfg$geneSets), 3), collapse = ", ")
    shiny::conditionalPanel("!output.hasResults", class = "nfr-welcome",
        bslib::card(class = "border-0 shadow-none bg-transparent",
            bslib::card_body(class = "text-center py-5",
                htmltools::img(src = "nfr-assets/logo.svg", height = 72,
                    alt = ""),
                htmltools::h2(class = "mt-3 fw-semibold",
                    "Which pathways is your gene list close to?"),
                htmltools::p(class = "lead text-body-secondary mx-auto",
                    style = "max-width: 46rem;",
                    "NetFunrankR tests gene lists against gene sets in a",
                    "molecular interaction network and tells you whether the",
                    "evidence comes from annotated members, from genes that",
                    "interact with the gene set, or from both, with exact,",
                    "calibrated p-values."),
                bslib::layout_column_wrap(width = "15rem", class = "my-4",
                    fill = FALSE, .step(1, "list-ul", "Paste genes",
                        "From GWAS, differential expression or curation."),
                    .step(2, "collection", "Pick a library",
                        paste0(libs, ", or your own GMT file.")),
                    .step(3, "bar-chart-line", "Explore",
                        "Evidence map, linker networks, report.")),
                if (.hasExample(cfg)) shiny::actionButton("tryExample",
                    paste("Try the", cfg$exampleLabel, "example"),
                    icon = bsicons::bs_icon("play-circle"),
                    class = "btn-primary btn-lg"))))
}

.step <- function(i, icon, title, text) {
    bslib::card(class = "h-100 nfr-step", bslib::card_body(
        htmltools::div(class = "nfr-step-num", i),
        bsicons::bs_icon(icon, size = "1.6rem", class = "text-primary"),
        htmltools::h5(class = "mt-2 mb-1", title),
        htmltools::p(class = "small text-body-secondary mb-0", text)))
}

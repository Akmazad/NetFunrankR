.resultsUI <- function() {
    shiny::conditionalPanel("output.hasResults", class = "nfr-results",
        bslib::layout_column_wrap(width = "13rem", fill = FALSE,
            class = "mb-3",
            .valueBox("Query genes in network", "vbGenes", "people",
                "primary"),
            .valueBox("Gene sets tested", "vbSets", "collection"),
            .valueBox("Significant gene sets", "vbSig", "check2-circle",
                "success"),
            .valueBox("Supported by network only", "vbNet", "diagram-3",
                "warning")),
        bslib::navset_card_underline(id = "tabs", full_screen = TRUE,
            bslib::nav_panel("Evidence map", value = "map",
                .help(paste("Each point is a gene set. x: overlap with the",
                    "query (degree-aware); y: linker genes beyond the",
                    "overlap. Click a point to open it in the explorer.")),
                plotly::plotlyOutput("evidenceMap", height = "520px")),
            bslib::nav_panel("Top gene sets", value = "top",
                plotly::plotlyOutput("topSets", height = "560px")),
            bslib::nav_panel("Results table", value = "table",
                reactable::reactableOutput("table")),
            bslib::nav_panel("Pathway explorer", value = "explorer",
                .explorerUI()),
            bslib::nav_panel("Compare with ORA", value = "ora",
                .help(paste("Hypergeometric ORA versus NetFunrankR for every",
                    "gene set. Points above the diagonal gain from network",
                    "evidence.")),
                plotly::plotlyOutput("oraCompare", height = "500px")),
            bslib::nav_panel("Reproduce in R", value = "code",
                .codeUI()),
            bslib::nav_spacer(),
            bslib::nav_item(.downloadMenu())))
}

.valueBox <- function(title, id, icon, theme = NULL) {
    bslib::value_box(title = title, value = shiny::textOutput(id),
        showcase = bsicons::bs_icon(icon, size = "1.5rem"),
        showcase_layout = "top right",
        theme = if (!is.null(theme)) paste0("text-", theme),
        class = "nfr-vbox border", shiny::uiOutput(paste0(id, "Sub")))
}

.help <- function(text) {
    htmltools::p(class = "small text-body-secondary mb-2",
        bsicons::bs_icon("info-circle", class = "me-1"), text)
}

.explorerUI <- function() {
    bslib::layout_columns(
        col_widths = bslib::breakpoints(sm = 12, lg = c(4, 8)),
        bslib::card(
            bslib::card_header("Gene set"),
            shiny::selectizeInput("pathway", NULL, choices = NULL,
                width = "100%", options = list(placeholder = "Search\u2026")),
            shiny::uiOutput("pathwayStats")),
        bslib::card(full_screen = TRUE,
            bslib::card_header(class = "d-flex flex-wrap gap-2",
                "Query genes linked to the gene set",
                htmltools::span(class = "ms-auto d-flex gap-2",
                    .legendPill("Query gene in set", "overlap + network"),
                    .legendPill("Linker gene", "network only"),
                    htmltools::span(class = "nfr-pill",
                        style = "--pill: #9fb6c3;", "Set member"))),
            visNetwork::visNetworkOutput("linkerNet", height = "460px"),
            bslib::card_footer(
                plotly::plotlyOutput("nullPlot", height = "190px"))))
}

.legendPill <- function(label, evidence) {
    htmltools::span(class = "nfr-pill",
        style = sprintf("--pill: %s;", .evidenceCols[[evidence]]), label)
}

.codeUI <- function() {
    htmltools::tagList(
        .help("Run this code in R to reproduce the analysis exactly."),
        htmltools::div(class = "position-relative",
            htmltools::tags$button(type = "button", id = "copyCode",
                class = "btn btn-sm btn-outline-secondary nfr-copy",
                `data-target` = "rcode", bsicons::bs_icon("clipboard"),
                "Copy"),
            shiny::verbatimTextOutput("rcode")))
}

.downloadMenu <- function() {
    bslib::popover(placement = "bottom",
        htmltools::tags$button(type = "button",
            class = "btn btn-sm btn-outline-primary my-1",
            bsicons::bs_icon("download"), "Download"),
        title = "Download",
        htmltools::div(class = "d-grid gap-2",
            shiny::downloadButton("dlCsv", "Results table (CSV)",
                class = "btn-sm"),
            shiny::downloadButton("dlReport", "HTML report", class = "btn-sm"),
            shiny::downloadButton("dlCode", "R script", class = "btn-sm")))
}

.guideUI <- function() {
    bslib::layout_columns(col_widths = bslib::breakpoints(sm = 12,
        lg = c(-2, 8, -2)),
        bslib::card(bslib::card_body(class = "nfr-prose",
            htmltools::h3("How NetFunrankR works"),
            htmltools::p("Over-representation analysis only counts query",
                "genes that are annotated to a pathway. NetFunrankR also",
                "asks whether the remaining query genes sit next to the",
                "pathway in the protein interaction network."),
            htmltools::tags$ol(
                htmltools::tags$li(htmltools::strong("Overlap test."),
                    "How many query genes are pathway members, compared",
                    "with random genes of the same network degree?"),
                htmltools::tags$li(htmltools::strong("Linker test."),
                    "Among the query genes that are not members, how many",
                    "interact directly with a member?"),
                htmltools::tags$li(htmltools::strong("Combined test."),
                    "The two tests are independent under the null",
                    "hypothesis and are combined with Fisher's method.")),
            htmltools::p("All p-values are computed exactly from the degree",
                "distribution of the network, without permutations."),
            htmltools::h4("Reading the results"),
            htmltools::tags$ul(
                .evidenceItem("overlap + network", "both tests significant"),
                .evidenceItem("overlap only", "driven by annotated members"),
                .evidenceItem("network only", paste("driven by linker genes:",
                    "associations that annotation-based enrichment misses")),
                .evidenceItem("joint", paste("neither test alone, only the",
                    "combination, is significant"))),
            htmltools::h4("Tips"),
            htmltools::tags$ul(
                htmltools::tags$li("Use HGNC symbols; unmatched genes are",
                    "listed under the gene box."),
                htmltools::tags$li("For differential-expression lists,",
                    "upload all measured genes as background."),
                htmltools::tags$li("Lists of a few genes work: network",
                    "evidence helps most when overlap is sparse.")))))
}

.evidenceItem <- function(label, text) {
    htmltools::tags$li(htmltools::span(class = "nfr-pill",
        style = sprintf("--pill: %s;", .evidenceCols[[label]]), label),
        " ", text)
}

.aboutUI <- function(cfg = .getAppConfig()) {
    cite <- tryCatch(format(utils::citation("NetFunrankR"), style = "text"),
        error = function(e) "See citation(\"NetFunrankR\").")
    bslib::layout_columns(col_widths = bslib::breakpoints(sm = 12,
        lg = c(-2, 8, -2)),
        bslib::card(bslib::card_body(class = "nfr-prose",
            htmltools::h3("About"),
            htmltools::p("NetFunrankR is an R/Bioconductor package for",
                "decomposed network-proximity enrichment analysis. This",
                "dashboard runs the same functions as the package; the",
                "'Reproduce in R' tab gives the code for every analysis."),
            htmltools::p(htmltools::strong("Citation: "), cite[1]),
            htmltools::p(htmltools::strong("Gene-set libraries: "),
                paste(names(cfg$geneSets), collapse = "; "), "."),
            htmltools::p(htmltools::strong("Interactomes: "),
                paste(names(cfg$networks), collapse = "; "), "."),
            htmltools::p(htmltools::strong("Version: "),
                as.character(utils::packageVersion("NetFunrankR"))),
            htmltools::p(htmltools::strong("Privacy: "), "uploaded gene lists",
                "and files are processed in memory for your session and are",
                "not stored by the application."))))
}

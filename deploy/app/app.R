# Entry point of the NetFunrankR dashboard for hosted deployments
# (Hugging Face Spaces / Docker, shinyapps.io, Posit Connect Cloud, Shiny
# Server). Configure with environment variables instead of editing code:
#
#   NFR_WORKERS        background worker processes        (default 1)
#   NFR_ASYNC          "false" computes in the app process instead of
#                      background workers (lower memory on small servers)
#   NFR_MAX_UPLOAD_MB  maximum upload size in MB          (default 100)
#   NFR_LIBRARIES      comma-separated subset of bundled libraries to offer,
#                      e.g. "KEGG 2021,WikiPathways 2024" (default: all)
#   NFR_GMT_DIR        folder of extra .gmt files offered as libraries
#   NFR_EDGE_DIR       folder of extra edge lists (.csv/.tsv) offered as
#                      networks
library(NetFunrankR)

env <- function(name, default) {
    v <- Sys.getenv(name, "")
    if (nzchar(v)) v else default
}

defaults <- dashboardDefaults()
libs <- defaults$geneSets
pick <- trimws(strsplit(env("NFR_LIBRARIES", ""), ",")[[1]])
if (length(pick)) libs <- libs[intersect(pick, names(libs))]
gmtDir <- env("NFR_GMT_DIR", "")
if (nzchar(gmtDir)) {
    gmt <- list.files(gmtDir, pattern = "\\.gmt$", full.names = TRUE)
    libs <- c(libs, stats::setNames(as.list(gmt),
        tools::file_path_sans_ext(basename(gmt))))
}

nets <- defaults$networks
edgeDir <- env("NFR_EDGE_DIR", "")
if (nzchar(edgeDir)) {
    ed <- list.files(edgeDir, pattern = "\\.(csv|tsv)$", full.names = TRUE)
    nets <- c(nets, stats::setNames(lapply(ed, function(f) {
        force(f)
        function() buildNetwork(utils::read.delim(f,
            sep = if (grepl("\\.csv$", f)) "," else "\t"), source = basename(f))
    }), tools::file_path_sans_ext(basename(ed))))
}

options(NetFunrankR.async = tolower(env("NFR_ASYNC", "true")) == "true")

netFunrankRApp(geneSets = libs, networks = nets,
    workers = as.integer(env("NFR_WORKERS", "1")),
    maxUploadMB = as.numeric(env("NFR_MAX_UPLOAD_MB", "100")))

# One-off migration of method labels in saved results to the final naming
# used in the manuscript (no recomputation; labels only):
#   NFR-combined       -> NetFunrankR
#   NFR-<component>    -> NetFunrankR-<component>
#   NetFunrankR-v1     -> Proximity-symmetric   (symmetric module distance)
#   Proximity-z        -> Proximity-closest     (closest distance, Guney 2016)
suppressPackageStartupMessages(library(data.table))
RES <- file.path(Sys.getenv("NF_ANALYSIS", "analysis"), "results")

## exact replacements, longest keys first
map <- c(
  "NFR-combined-dist" = "NetFunrankR-distance-combined",
  "NFR-combined" = "NetFunrankR",
  "NFR-overlap" = "NetFunrankR-overlap",
  "NFR-linker" = "NetFunrankR-linker",
  "NFR-distance" = "NetFunrankR-distance",
  "NFR-beyond" = "NetFunrankR-distance",
  "NFR-total" = "NetFunrankR-total",
  "NetFunrankR-v1" = "Proximity-symmetric",
  "Proximity-z" = "Proximity-closest",
  "NetFunrankR v1 (original code)" =
    "Original implementation (symmetric proximity, 100 perm)",
  "NFR (all four tests)" = "NetFunrankR (all tests)",
  "Proximity-z (1000 perm)" = "Proximity-closest (1000 perm)",
  "NetFunrankR-v1 statistic (1000 perm)" = "Proximity-symmetric (1000 perm)",
  "NFR" = "NetFunrankR")

relabel <- function(v) {
  hit <- !is.na(v) & v %in% names(map)
  v[hit] <- map[v[hit]]
  v
}

fix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    for (j in seq_along(x)) {
      if (is.character(x[[j]])) x[[j]] <- relabel(x[[j]])
      if (is.factor(x[[j]])) levels(x[[j]]) <- relabel(levels(x[[j]]))
    }
    names(x) <- relabel(names(x))
    return(x)
  }
  if (is.matrix(x)) {
    if (!is.null(colnames(x))) colnames(x) <- relabel(colnames(x))
    return(x)
  }
  if (is.list(x) && !methods::is(x, "DataFrame")) {
    if (!is.null(names(x))) names(x) <- relabel(names(x))
    return(lapply(x, fix))
  }
  if (is.character(x) && !is.null(names(x))) names(x) <- relabel(names(x))
  if (is.numeric(x) && !is.null(names(x))) names(x) <- relabel(names(x))
  x
}

files <- list.files(RES, pattern = "\\.(rds|csv)$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("proto_adjacent", files)]
for (f in files) {
  if (grepl("\\.csv$", f)) {
    d <- fread(f)
    d2 <- as.data.table(fix(d))
    if (!identical(d, d2)) {
      fwrite(d2, f)
      cat("relabelled", f, "\n")
    }
  } else {
    o <- readRDS(f)
    keepDT <- is.data.table(o)
    o2 <- fix(o)
    if (keepDT) o2 <- as.data.table(o2)
    if (!identical(o, o2)) {
      saveRDS(o2, f)
      cat("relabelled", f, "\n")
    }
  }
}
unlink(file.path(RES, "gwas", "proto_adjacent.rds"))

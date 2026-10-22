# Build interactome variants from STRING v12 per-channel evidence.
#
#   string_full700   : STRING combined score >= 700 (all channels; contains the
#                      curated-database and text-mining channels -> circular
#                      with pathway annotations, used only as sensitivity)
#   string_nokb700   : combined score recomputed WITHOUT the "database" and
#                      "textmining" channels, >= 700  (primary network)
#   string_nokb400   : same, >= 400 (sensitivity)
#
# Output: analysis/data/networks.rds  (named list of 2-column symbol edge lists)

suppressPackageStartupMessages({
  library(data.table)
})
root <- Sys.getenv("NF_ANALYSIS", "analysis")
raw  <- file.path(root, "data", "raw")

info <- fread(file.path(raw, "9606.protein.info.v12.0.txt.gz"),
              select = 1:2, col.names = c("id", "symbol"))
lk <- fread(file.path(raw, "9606.protein.links.detailed.v12.0.txt.gz"))
cat("STRING detailed links:", nrow(lk), "\n")

# STRING channel-combination (von Mering et al. 2005; Szklarczyk et al. 2021):
# remove the prior from each channel, combine as independent evidence,
# then add the prior back.
string_combine <- function(dt, channels, prior = 0.041) {
  acc <- rep(1, nrow(dt))
  for (ch in channels) {
    s <- dt[[ch]] / 1000
    s <- pmax((s - prior) / (1 - prior), 0)
    acc <- acc * (1 - s)
  }
  comb <- 1 - acc
  round(1000 * (comb + prior * (1 - comb)))
}

nokb_channels <- c("neighborhood", "fusion", "cooccurence",
                   "coexpression", "experimental")
lk[, nokb_score := string_combine(.SD, nokb_channels), .SDcols = nokb_channels]

lk[, from := info$symbol[match(protein1, info$id)]]
lk[, to   := info$symbol[match(protein2, info$id)]]
lk <- lk[!is.na(from) & !is.na(to) & from != to]
# undirected: keep each pair once
lk <- lk[from < to]

mk <- function(score, th) unique(lk[get(score) >= th, .(from, to)])
nets <- list(
  string_full700 = mk("combined_score", 700),
  string_nokb700 = mk("nokb_score", 700),
  string_nokb400 = mk("nokb_score", 400)
)
for (n in names(nets)) {
  e <- nets[[n]]
  cat(sprintf("%-16s edges=%7d nodes=%6d\n", n, nrow(e),
              length(unique(c(e$from, e$to)))))
}
saveRDS(nets, file.path(root, "data", "networks.rds"))

# Build interactome variants from STRING per-channel evidence.
#
#   string_nokb400 : combined score recomputed WITHOUT the "database" and
#                    "textmining" channels, >= 400  (primary network)
#   string_nokb700 : same, >= 700 (sensitivity)
#   string_full700 : published combined score >= 700 (all channels; contains
#                    curated-database and text-mining evidence, i.e. circular
#                    with pathway annotations; sensitivity only)
#
# Settings (environment variables):
#   NF_STRING_VERSION  STRING version            (default "12.0")
#   NF_TAXON           NCBI taxonomy identifier  (default "9606", human)
# Missing source files are downloaded to analysis/data/raw.
# Output: analysis/data/networks.rds  (named list of 2-column edge lists)

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("NF_ANALYSIS", "analysis")
raw <- file.path(root, "data", "raw")
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
version <- Sys.getenv("NF_STRING_VERSION", "12.0")
taxon <- Sys.getenv("NF_TAXON", "9606")

string_file <- function(kind) {
  f <- sprintf("%s.%s.v%s.txt.gz", taxon, kind, version)
  path <- file.path(raw, f)
  if (!file.exists(path)) {
    url <- sprintf("https://stringdb-downloads.org/download/%s.v%s/%s",
                   kind, version, f)
    message("downloading ", url)
    options(timeout = max(3600, getOption("timeout")))
    utils::download.file(url, path, mode = "wb")
  }
  path
}

info <- fread(string_file("protein.info"), select = 1:2,
              col.names = c("id", "symbol"))
lk <- fread(string_file("protein.links.detailed"))
cat("STRING detailed links:", nrow(lk), "\n")

# STRING channel-combination rule (von Mering et al. 2005): remove the prior
# from each channel, combine as independent evidence, add the prior back.
string_combine <- function(dt, channels, prior = 0.041) {
  acc <- rep(1, nrow(dt))
  for (ch in channels) {
    s <- pmax((dt[[ch]] / 1000 - prior) / (1 - prior), 0)
    acc <- acc * (1 - s)
  }
  comb <- 1 - acc
  round(1000 * (comb + prior * (1 - comb)))
}

nokb_channels <- c("neighborhood", "fusion", "cooccurence",
                   "coexpression", "experimental")
lk[, nokb_score := string_combine(.SD, nokb_channels), .SDcols = nokb_channels]
lk[, from := info$symbol[match(protein1, info$id)]]
lk[, to := info$symbol[match(protein2, info$id)]]
lk <- lk[!is.na(from) & !is.na(to) & from != to]
lk <- lk[from < to]   # undirected: keep each pair once

variants <- data.table(name = c("string_nokb400", "string_nokb700", "string_full700"),
                       score = c("nokb_score", "nokb_score", "combined_score"),
                       threshold = c(400, 700, 700))
nets <- setNames(lapply(seq_len(nrow(variants)), function(i)
  unique(lk[get(variants$score[i]) >= variants$threshold[i], .(from, to)])),
  variants$name)
for (n in names(nets))
  cat(sprintf("%-16s edges=%7d nodes=%6d\n", n, nrow(nets[[n]]),
              length(unique(c(nets[[n]]$from, nets[[n]]$to)))))
saveRDS(nets, file.path(root, "data", "networks.rds"))

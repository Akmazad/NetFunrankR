# Provenance of data/STRING_v12_PPI.rda
#
# Source files (STRING v12.0, https://string-db.org, CC BY 4.0):
#   https://stringdb-downloads.org/download/protein.links.detailed.v12.0/
#       9606.protein.links.detailed.v12.0.txt.gz
#   https://stringdb-downloads.org/download/protein.info.v12.0/
#       9606.protein.info.v12.0.txt.gz
#
# The STRING combined score is recomputed without the "database" and
# "textmining" channels using STRING's channel-combination rule
# (prior 0.041; von Mering et al. 2005). Proteins are mapped to their STRING
# preferred names (HGNC symbols). Edges are kept if the recomputed score is
# >= 400 or the published combined score is >= 700.

library(data.table)
raw <- "raw"   # directory with the two downloaded files
info <- fread(file.path(raw, "9606.protein.info.v12.0.txt.gz"),
    select = 1:2, col.names = c("id", "symbol"))
lk <- fread(file.path(raw, "9606.protein.links.detailed.v12.0.txt.gz"))

combine <- function(dt, channels, prior = 0.041) {
    acc <- rep(1, nrow(dt))
    for (ch in channels) {
        s <- pmax((dt[[ch]] / 1000 - prior) / (1 - prior), 0)
        acc <- acc * (1 - s)
    }
    comb <- 1 - acc
    as.integer(round(1000 * (comb + prior * (1 - comb))))
}
channels <- c("neighborhood", "fusion", "cooccurence", "coexpression",
    "experimental")
lk[, nokb_score := combine(.SD, channels), .SDcols = channels]
lk[, from := info$symbol[match(protein1, info$id)]]
lk[, to := info$symbol[match(protein2, info$id)]]
lk <- lk[!is.na(from) & !is.na(to) & from < to &
    (nokb_score >= 400 | combined_score >= 700)]
lk <- lk[, .(combined_score = max(combined_score),
    nokb_score = max(nokb_score)), by = .(from, to)]
STRING_v12_PPI <- as.data.frame(lk)
save(STRING_v12_PPI, file = "STRING_v12_PPI.rda", compress = "xz")

# Pathway libraries (KEGG_2021_Human, Reactome_Pathways_2024,
# WikiPathways_2024_Human) were downloaded from the Enrichr library page
# (https://maayanlab.cloud/Enrichr/#libraries) and converted to two-column
# data.frames (GeneSymbol, PathwayName). COPD_proteinCodingGenes are the
# protein-coding COPD genes of Morselli Gysi et al. (2021, PNAS 118:
# e2025581118, Dataset S1).

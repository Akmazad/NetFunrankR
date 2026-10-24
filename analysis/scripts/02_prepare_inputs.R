# Prepare all benchmark inputs (run after 01_build_networks.R).
#
#   data/nf_networks.rds      NFNetwork objects for the network variants
#   data/kegg_gs_entrez.rds   KEGG gene sets (EnrichmentBrowser::getGenesets)
#   data/geo2kegg_de.rds      GEO2KEGG compendium, summarised + limma DE
#   data/gwas_sig.rds         genome-wide significant GWAS Catalog hits
#
# Settings (environment variables):
#   NF_GWAS_URL  GWAS Catalog associations archive (default: latest release)
#   NF_GWAS_P    significance threshold for GWAS hits (default 5e-8)
source("analysis/scripts/00_setup.R")
dat <- file.path(AN, "data")
raw <- file.path(dat, "raw")
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(3600, getOption("timeout")))

# ---------------------------------------------------------------- networks
saveRDS(load_networks(), file.path(dat, "nf_networks.rds"))

# ---------------------------------------------------------------- KEGG gene sets
f_gs <- file.path(dat, "kegg_gs_entrez.rds")
if (!file.exists(f_gs))
  saveRDS(EnrichmentBrowser::getGenesets(org = "hsa", db = "kegg"), f_gs)

# ---------------------------------------------------------------- GEO2KEGG
f_exp <- file.path(dat, "geo2kegg_de.rds")
if (!file.exists(f_exp)) {
  suppressPackageStartupMessages(library(GSEABenchmarkeR))
  geo2kegg <- loadEData("geo2kegg")
  geo2kegg <- maPreproc(geo2kegg)
  geo2kegg <- runDE(geo2kegg, de.method = "limma", padj.method = "flexible")
  saveRDS(geo2kegg, f_exp)
}

# ---------------------------------------------------------------- GWAS Catalog
gwas_url <- Sys.getenv("NF_GWAS_URL", paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/",
  "gwas-catalog-associations_ontology-annotated-full.zip"))
gwas_p <- as.numeric(Sys.getenv("NF_GWAS_P", "5e-8"))
f_gwas <- file.path(dat, "gwas_sig.rds")
if (!file.exists(f_gwas)) {
  zip <- file.path(raw, basename(gwas_url))
  if (!file.exists(zip)) utils::download.file(gwas_url, zip, mode = "wb")
  tsv <- list.files(raw, pattern = "^gwas-catalog-download-associations.*\\.tsv$",
                    full.names = TRUE)
  if (!length(tsv)) {
    utils::unzip(zip, exdir = raw)
    tsv <- list.files(raw, pattern = "^gwas-catalog-download-associations.*\\.tsv$",
                      full.names = TRUE)
  }
  g <- fread(tsv[1], quote = "",
             select = c("DATE", "MAPPED_TRAIT", "MAPPED_GENE", "PVALUE_MLOG", "PUBMEDID"))
  g <- g[PVALUE_MLOG >= -log10(gwas_p)]
  attr(g, "source") <- list(file = basename(tsv[1]), url = gwas_url,
                            p_threshold = gwas_p, retrieved = Sys.Date())
  saveRDS(g, f_gwas)
  cat("GWAS hits at p <=", gwas_p, ":", nrow(g), "\n")
}

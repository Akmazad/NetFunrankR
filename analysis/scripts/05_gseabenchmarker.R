# GSEABenchmarkeR benchmark (Geistlinger et al. 2021, Brief Bioinform):
# 42 GEO microarray datasets (geo2kegg), KEGG gene sets, MalaCards relevance.
#
# Comparators (EnrichmentBrowser implementations):
#   ora (hypergeometric), camera, gsva, globaltest        [perm = 0]
#   safe, gsea, padog, roast, samgs                         [perm = 1000]
# Gene-list network methods (on the DE genes used by ORA, ADJ.PVAL < 0.05):
#   NEA, RWR, Proximity-closest (Guney et al. 2016), NetFunrankR (this work)
#
# Output: analysis/results/gseabench/  (rankings, runtimes, relevance scores)
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GSEABenchmarkeR)
  library(EnrichmentBrowser)
  library(BiocParallel)
  library(org.Hs.eg.db)
})
OUT <- file.path(RES, "gseabench")
dir.create(OUT, showWarnings = FALSE)
cores <- nf_cores()

# ---------------------------------------------------------------- data
f_exp <- file.path(AN, "data", "geo2kegg_de.rds")
if (!file.exists(f_exp)) {
  geo2kegg <- loadEData("geo2kegg")
  geo2kegg <- maPreproc(geo2kegg)
  geo2kegg <- runDE(geo2kegg, de.method = "limma", padj.method = "flexible")
  saveRDS(geo2kegg, f_exp)
}
geo2kegg <- readRDS(f_exp)
ds <- Sys.getenv("NF_DATASETS", "")
if (nzchar(ds)) geo2kegg <- geo2kegg[strsplit(ds, ",")[[1]]]

f_gs <- file.path(AN, "data", "kegg_gs_entrez.rds")
if (!file.exists(f_gs)) saveRDS(getGenesets(org = "hsa", db = "kegg"), f_gs)
kegg.gs <- readRDS(f_gs)

# ---------------------------------------------------------------- EnrichmentBrowser methods
eb_fast <- c("ora", "camera", "gsva", "globaltest")
eb_perm <- c("safe", "gsea", "padog", "roast", "samgs")
run_eb <- Sys.getenv("NF_RUN_EB", "TRUE") == "TRUE"
if (run_eb) {
  bp <- SnowParam(workers = cores, type = "SOCK")
  runEA(geo2kegg, eb_fast, kegg.gs, perm = 0, parallel = bp,
        save2file = TRUE, out.dir = OUT)
  runEA(geo2kegg, eb_perm, kegg.gs, perm = 1000, parallel = bp,
        save2file = TRUE, out.dir = OUT)
}

# ---------------------------------------------------------------- network methods
eg2sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(c(unlist(kegg.gs),
                                  unlist(lapply(geo2kegg, rownames)))),
                                keytype = "ENTREZID", column = "SYMBOL")
nets <- readRDS(file.path(AN, "data", "nf_networks.rds"))
net <- nets$string_nokb400

# one task per dataset; gene sets restricted to measured genes, 5-500 genes
# (as in EnrichmentBrowser), mapped to symbols
tasks <- lapply(names(geo2kegg), function(id) {
  se <- geo2kegg[[id]]
  rd <- SummarizedExperiment::rowData(se)
  measured <- rownames(se)
  gs <- lapply(kegg.gs, intersect, measured)
  gs <- gs[lengths(gs) >= 5 & lengths(gs) <= 500]
  list(id = id,
       query = unname(stats::na.omit(eg2sym[measured[rd$ADJ.PVAL < 0.05]])),
       background = unname(stats::na.omit(eg2sym[measured])),
       gs = lapply(gs, function(s) unname(stats::na.omit(eg2sym[s]))))
})

cl <- parallel::makeCluster(cores)
parallel::clusterExport(cl, "net")
parallel::clusterEvalQ(cl, source("analysis/scripts/00_setup.R"))
parallel::clusterSetRNGStream(cl, 42)
net_res <- parallel::parLapplyLB(cl, tasks, function(t) {
  ctx <- make_ctx(t$gs, net, min_size = 5L, max_size = 500L)
  q <- t$query
  out <- list()
  tm <- c()
  nfr_names <- c("NetFunrankR-overlap", "NetFunrankR-linker", "NetFunrankR-distance", "NetFunrankR-total",
                 "NetFunrankR", "NetFunrankR-distance-combined")
  if (sum(nodeNames(net) %in% q) < 2L) {
    p1 <- stats::setNames(rep(1, length(ctx$keep)), ctx$keep)
    for (m in c("NEA", "RWR", "Proximity-closest", nfr_names)) out[[m]] <- p1
    return(list(id = t$id, p = out, time = tm, n_query = length(q)))
  }
  tm["NEA"] <- system.time(out[["NEA"]] <- m_nea(q, ctx))[3]
  tm["RWR"] <- system.time(out[["RWR"]] <- m_rwr(q, ctx))[3]
  tm["Proximity-closest"] <- system.time(out[["Proximity-closest"]] <- m_prox_closest(q, ctx))[3]
  # timing of one NetFunrankR run (m_nfr runs the test twice: linker + distance)
  tm["NetFunrankR"] <- system.time(r <- m_nfr(q, ctx, background = t$background))[3] / 2
  out[nfr_names] <- r[nfr_names]
  list(id = t$id, p = out, time = tm, n_query = length(q))
})
parallel::stopCluster(cl)
saveRDS(net_res, file.path(OUT, "network_methods_raw.rds"))

# ---------------------------------------------------------------- rankings
# As in EnrichmentBrowser rankings, only tested sets are ranked (untested
# sets contribute nothing to the relevance score).
to_ranking <- function(p, all_sets) {
  p <- p[!is.na(p)]
  p <- p[order(p)]
  S4Vectors::DataFrame(GENE.SET = names(p), PVAL = unname(p))
}
eb_methods <- c(eb_fast, eb_perm)
ea.ranks <- readResults(OUT, names(geo2kegg), methods = eb_methods, type = "ranking")
ea.rtimes <- readResults(OUT, names(geo2kegg), methods = eb_methods, type = "runtime")
all_sets <- names(kegg.gs)
for (m in names(net_res[[1]]$p)) {
  ea.ranks[[m]] <- stats::setNames(lapply(net_res, function(x) to_ranking(x$p[[m]], all_sets)),
                                   vapply(net_res, `[[`, "", "id"))
}
saveRDS(ea.ranks, file.path(OUT, "all_rankings.rds"))
saveRDS(list(eb = ea.rtimes, net = lapply(net_res, `[[`, "time")),
        file.path(OUT, "runtimes.rds"))

# ---------------------------------------------------------------- evaluation
mala.kegg.file <- system.file("extdata/malacards/KEGG.rds", package = "GSEABenchmarkeR")
mala.kegg <- readRDS(mala.kegg.file)
d2d.file <- system.file("extdata/malacards/GseId2Disease.txt", package = "GSEABenchmarkeR")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)
d2d.map <- d2d.map[names(geo2kegg)]

rel <- evalRelevance(ea.ranks, mala.kegg, d2d.map)     # % of optimal score
saveRDS(rel, file.path(OUT, "relevance.rds"))
nsig <- evalNrSigSets(ea.ranks, alpha = 0.05, padj = "BH")
saveRDS(nsig, file.path(OUT, "nr_sig_sets.rds"))
print(round(sort(apply(rel, 2, mean), decreasing = TRUE), 2))
print(round(sort(apply(nsig, 2, median), decreasing = TRUE), 2))

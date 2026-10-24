# Case study: COPD genes (Morselli Gysi et al. 2021, PNAS) vs KEGG.
#  - full NetFunrankR result table + evidence classes
#  - agreement with MalaCards COPD (PDCO) relevance, all gene-list methods
#  - cross-modal check: CAMERA on COPD expression data (GSE42057, geo2kegg)
#  - small-list scenario: random subsets of 6 / 10 / 20 COPD genes
# Output: analysis/results/copd/
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GSEABenchmarkeR)
  library(parallel)
})
OUT <- file.path(RES, "copd")
dir.create(OUT, showWarnings = FALSE)
net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
q <- COPD_proteinCodingGenes

kegg.gs <- readRDS(file.path(AN, "data", "kegg_gs_entrez.rds"))
eg2sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(unlist(kegg.gs)),
                                "SYMBOL", "ENTREZID")
gs <- lapply(kegg.gs, function(s) unique(unname(stats::na.omit(eg2sym[s]))))
ctx <- make_ctx(gs, net)

resObj <- netEnrich(q, gs, net, distances = ctx$D, maxSize = 500L)
saveRDS(resObj, file.path(OUT, "copd_kegg_result.rds"))
res <- as.data.frame(resObj)
fwrite(res, file.path(OUT, "copd_kegg_netfunrankr.csv"))
print(table(res$evidence))

set.seed(2)
pv <- run_genelist_methods(q, ctx)
saveRDS(pv, file.path(OUT, "copd_all_methods.rds"))

mala <- readRDS(system.file("extdata/malacards/KEGG.rds", package = "GSEABenchmarkeR"))$PDCO
opt <- compOpt(mala, unique(sub("_.*", "", pv$set)))  # compOpt expects bare KEGG ids
rel <- pv[!is.na(p), {
  dd <- .SD[order(p)]
  .(relevance = 100 * evalRelevance(S4Vectors::DataFrame(GENE.SET = dd$set, PVAL = dd$p), mala) / opt,
    n_fdr = sum(p.adjust(p, "BH") < 0.05))
}, by = method]
fwrite(rel, file.path(OUT, "copd_malacards_relevance.csv"))
print(rel[order(-relevance)])

# cross-modal: CAMERA on GSE42057 (COPD PBMC expression)
geo <- readRDS(file.path(AN, "data", "geo2kegg_de.rds"))[["GSE42057"]]
cam <- EnrichmentBrowser::sbea("camera", geo, kegg.gs, perm = 0)
cam <- as.data.frame(EnrichmentBrowser::gsRanking(cam, signif.only = FALSE))
cam_p <- stats::setNames(cam$PVAL, cam$GENE.SET)
xm <- pv[!is.na(p), .(rho = suppressWarnings(stats::cor(-log10(p), -log10(cam_p[set]),
                                                       method = "spearman", use = "complete.obs"))),
         by = method]
fwrite(xm, file.path(OUT, "copd_vs_GSE42057_camera.csv"))
print(xm[order(-rho)])

# small query lists: random subsets of COPD genes
cl <- makeCluster(nf_cores())
clusterExport(cl, c("gs", "net", "q", "mala", "opt"))
clusterEvalQ(cl, {
  source("analysis/scripts/00_setup.R")
  ctx <- make_ctx(gs, net)
  NULL
})
clusterSetRNGStream(cl, 6)
jobs <- CJ(size = c(6L, 10L, 20L, 50L), rep = 1:50)
small <- rbindlist(parLapplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) {
  j <- jobs[i]
  qq <- sample(intersect(q, nodeNames(net)), j$size)
  r <- run_genelist_methods(qq, ctx, methods = c("ora", "nea", "rwr", "closest", "nfr"))
  r[!is.na(p), {
    dd <- .SD[order(p)]
    .(relevance = 100 * GSEABenchmarkeR::evalRelevance(
      S4Vectors::DataFrame(GENE.SET = dd$set, PVAL = dd$p), mala) / opt,
      n_fdr = sum(p.adjust(p, "BH") < 0.05))
  }, by = method][, `:=`(size = j$size, rep = j$rep)]
}, jobs = jobs))
stopCluster(cl)
fwrite(small, file.path(OUT, "copd_small_lists.csv"))
print(dcast(small[, .(rel = mean(relevance)), by = .(method, size)], method ~ size, value.var = "rel"))

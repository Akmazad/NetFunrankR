# Type-I error calibration under random (null) query gene lists.
#   uniform : genes drawn uniformly from the network
#   hub     : genes drawn with probability proportional to degree (mimics the
#             hub/study bias of real disease-gene and GWAS lists)
# Output: analysis/results/calibration_pvalues.rds
source("analysis/scripts/00_setup.R")
library(parallel)

sizes <- c(20L, 50L, 100L, 200L)
schemes <- c("uniform", "hub")
n_rep <- as.integer(Sys.getenv("NF_NREP", 200))
design <- CJ(size = sizes, scheme = schemes, rep = seq_len(n_rep))

METHODS <- nf_methods_env(c("ora", "nea", "rwr", "guney", "v1", "nfr"))
cl <- makeCluster(as.integer(Sys.getenv("NF_CORES", 14)))
clusterExport(cl, "METHODS")
on.exit(stopCluster(cl))
clusterEvalQ(cl, {
  source("analysis/scripts/00_setup.R")
  nets <- readRDS(file.path(AN, "data", "nf_networks.rds"))
  net <- nets$string_nokb400
  ctx <- make_ctx(kegg_symbols(), net)
  NULL
})
clusterSetRNGStream(cl, 20261021)

res <- parLapplyLB(cl, seq_len(nrow(design)), function(i, design) {
  d <- design[i]
  w <- if (d$scheme == "hub") net$degree else rep(1, length(net$genes))
  q <- sample(net$genes, d$size, prob = w)
  r <- run_genelist_methods(q, ctx, methods = METHODS)
  r[, `:=`(size = d$size, scheme = d$scheme, rep = d$rep)]
  r
}, design = design)
res <- rbindlist(res)
f_out <- file.path(RES, "calibration_pvalues.rds")
if (nzchar(Sys.getenv("NF_METHODS")) && file.exists(f_out)) res <- merge_results(readRDS(f_out), res)
saveRDS(res, f_out)

summ <- res[!is.na(p), .(fpr05 = mean(p < 0.05), fpr01 = mean(p < 0.01)),
            by = .(method, scheme, size)]
fdr <- res[!is.na(p), .(n_fdr = sum(p.adjust(p, "BH") < 0.05)),
           by = .(method, scheme, size, rep)][
  , .(mean_fdr_hits = mean(n_fdr), frac_any = mean(n_fdr > 0)),
  by = .(method, scheme, size)]
summ <- merge(summ, fdr, by = c("method", "scheme", "size"))
fwrite(summ, file.path(RES, "calibration_summary.csv"))
print(dcast(summ, method ~ scheme + size, value.var = "fpr05"))
print(dcast(summ, method ~ scheme + size, value.var = "frac_any"))

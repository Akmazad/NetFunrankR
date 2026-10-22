# Power under incomplete annotation (planted-signal simulation).
# For a target pathway Y, a fraction h of its genes is hidden from the
# annotation that is tested (Y' = Y \ H). The query contains k = s * n genes
# drawn from the complete pathway Y plus n - k random genes. Hidden signal
# genes are invisible to overlap-based tests and can only be detected through
# their interactome proximity to the remaining annotated members.
# Output: analysis/results/power_holdout.rds
source("analysis/scripts/00_setup.R")
library(parallel)

n_query <- 50L
signal <- c(0.1, 0.2, 0.4)
hidden <- c(0, 0.25, 0.5, 0.75, 0.9)
n_target <- as.integer(Sys.getenv("NF_NTARGET", 100))

nets <- readRDS(file.path(AN, "data", "nf_networks.rds"))
net <- nets$string_nokb400
gs <- kegg_symbols()
n_in <- vapply(gs, function(s) sum(net$genes %in% s), 1L)
set.seed(7)
targets <- sample(names(n_in)[n_in >= 50 & n_in <= 200], n_target)
design <- CJ(target = targets, s = signal, h = hidden)
design[, id := .I]

METHODS <- nf_methods_env(c("ora", "nea", "rwr", "guney", "v1", "nfr"))
cl <- makeCluster(as.integer(Sys.getenv("NF_CORES", 14)))
clusterExport(cl, "METHODS")
on.exit(stopCluster(cl))
clusterExport(cl, "n_query")
clusterEvalQ(cl, {
  source("analysis/scripts/00_setup.R")
  nets <- readRDS(file.path(AN, "data", "nf_networks.rds"))
  net <- nets$string_nokb400
  gs <- kegg_symbols()
  ctx0 <- make_ctx(gs, net)
  # replace one gene set in a context (distances + incidence)
  ctx_replace <- function(ctx, name, genes) {
    ctx$gs[[name]] <- genes
    ctx$D[, name] <- set_distances(ctx$net, ctx$gs[name])[, 1]
    j <- match(name, ctx$keep)
    if (!is.na(j)) ctx$M[, j] <- as.numeric(ctx$net$genes %in% genes)
    ctx
  }
  NULL
})
clusterSetRNGStream(cl, 1103)

res <- parLapplyLB(cl, seq_len(nrow(design)), function(i, design) {
  d <- design[i]
  Y <- intersect(gs[[d$target]], net$genes)
  H <- sample(Y, round(d$h * length(Y)))
  ctx <- ctx_replace(ctx0, d$target, setdiff(gs[[d$target]], H))
  k <- round(d$s * n_query)
  sig <- sample(Y, k)
  bg <- sample(setdiff(net$genes, Y), n_query - k)
  r <- run_genelist_methods(c(sig, bg), ctx, methods = METHODS)
  r[, `:=`(q = p.adjust(p, "BH"), rank = rank(p, ties.method = "average") / .N), by = method]
  out <- r[set == d$target]
  out[, `:=`(target = d$target, s = d$s, h = d$h, id = d$id,
             n_hidden_signal = sum(sig %in% H), n_sets = nrow(r) / uniqueN(r$method))]
  out
}, design = design)
res <- rbindlist(res)
f_out <- file.path(RES, "power_holdout.rds")
if (nzchar(Sys.getenv("NF_METHODS")) && file.exists(f_out)) res <- merge_results(readRDS(f_out), res)
saveRDS(res, f_out)
summ <- res[, .(power = mean(q < 0.05), median_rank = median(rank)), by = .(method, s, h)]
fwrite(summ, file.path(RES, "power_holdout_summary.csv"))
print(dcast(summ[s == 0.2], method ~ h, value.var = "power"))

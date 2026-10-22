# Robustness of the GWAS-list benchmark to interactome choices.
#   network variants : STRING noKB >= 400 (primary), noKB >= 700, full >= 700
#   edge dropout     : 10 / 25 / 50 % of edges removed at random (3 reps)
#   rewired          : degree-preserving rewiring (negative control; 3 reps)
#   degree bins      : minimum bin size 50 / 100 / 250
#   subsampling      : 50 % of each disease list, 20 reps; Jaccard of top-20
# Output: analysis/results/robustness/
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GSEABenchmarkeR)
  library(parallel)
})
OUT <- file.path(RES, "robustness")
dir.create(OUT, showWarnings = FALSE)
raw <- readRDS(file.path(RES, "gwas", "gwas_raw.rds"))
lists <- lapply(raw$lists, `[[`, "genes")
kegg.gs <- readRDS(file.path(AN, "data", "kegg_gs_entrez.rds"))
eg2sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(unlist(kegg.gs)),
                                "SYMBOL", "ENTREZID")
gs <- lapply(kegg.gs, function(s) unique(unname(stats::na.omit(eg2sym[s]))))
edges <- readRDS(file.path(AN, "data", "networks.rds"))

variants <- list(list(name = "STRING noKB>=400", edges = "string_nokb400"),
                 list(name = "STRING noKB>=700", edges = "string_nokb700"),
                 list(name = "STRING all>=700 (circular)", edges = "string_full700"))
for (f in c(0.1, 0.25, 0.5)) for (r in 1:3)
  variants[[length(variants) + 1]] <- list(name = sprintf("dropout %d%%", 100 * f),
                                           edges = "string_nokb400", drop = f, rep = r)
for (r in 1:3)
  variants[[length(variants) + 1]] <- list(name = "rewired", edges = "string_nokb400",
                                           rewire = TRUE, rep = r)
for (b in c(50L, 250L))
  variants[[length(variants) + 1]] <- list(name = sprintf("bin size %d", b),
                                           edges = "string_nokb400", bin = b)

cl <- makeCluster(min(length(variants), as.integer(Sys.getenv("NF_CORES", 14))))
clusterExport(cl, c("edges", "gs", "lists"))
clusterEvalQ(cl, source("analysis/scripts/00_setup.R"))
clusterSetRNGStream(cl, 5)
res <- parLapplyLB(cl, variants, function(v) {
  e <- edges[[v$edges]]
  if (!is.null(v$drop)) e <- e[sample(nrow(e), round((1 - v$drop) * nrow(e))), ]
  if (isTRUE(v$rewire)) {
    g <- igraph::graph_from_data_frame(e, directed = FALSE)
    g <- igraph::rewire(g, igraph::keeping_degseq(niter = 10 * igraph::ecount(g)))
    e <- igraph::as_data_frame(g)[, 1:2]
  }
  net <- nf_network(e, min_bin_size = if (is.null(v$bin)) 100L else v$bin)
  ctx <- make_ctx(gs, net)
  out <- rbindlist(lapply(names(lists), function(d) {
    r <- run_genelist_methods(lists[[d]], ctx, methods = c("ora", "guney", "nfr"))
    r[, disease := d]
  }))
  out[, `:=`(variant = v$name, rep = if (is.null(v$rep)) 1L else v$rep)]
  out
})
pv <- rbindlist(res)
saveRDS(pv, file.path(OUT, "robustness_pvalues.rds"))

mala <- readRDS(system.file("extdata/malacards/KEGG.rds", package = "GSEABenchmarkeR"))
rel <- pv[, {
  dd <- .SD[!is.na(p)][order(p)]
  .(rel = evalRelevance(S4Vectors::DataFrame(GENE.SET = dd$set, PVAL = dd$p),
                        mala[[disease[1]]]) /
      compOpt(mala[[disease[1]]], unique(sub("_.*", "", dd$set))) * 100)
}, by = .(variant, rep, method, disease)]
fwrite(rel, file.path(OUT, "robustness_relevance.csv"))
print(dcast(rel[, .(rel = mean(rel)), by = .(variant, method)], variant ~ method, value.var = "rel"))

# ---------------------------------------------------------------- query subsampling stability
net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
clusterExport(cl, "net")
clusterEvalQ(cl, {
  ctx <- make_ctx(gs, net)
  NULL
})
jobs <- CJ(disease = names(lists), rep = 0:20)
stab <- parLapplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) {
  j <- jobs[i]
  q <- lists[[j$disease]]
  if (j$rep > 0) q <- sample(q, ceiling(length(q) / 2))
  r <- run_genelist_methods(q, ctx, methods = c("ora", "rwr", "guney", "nfr"))
  r[, `:=`(disease = j$disease, rep = j$rep)]
}, jobs = jobs)
stopCluster(cl)
stab <- rbindlist(stab)
top <- function(dt, k = 20) dt[!is.na(p)][order(p)][seq_len(min(k, .N)), set]
jac <- stab[rep > 0, {
  m0 <- .BY$method
  d0 <- .BY$disease
  ref <- top(stab[rep == 0 & method == m0 & disease == d0])
  .(jaccard = length(intersect(top(.SD), ref)) / length(union(top(.SD), ref)))
}, by = .(method, disease, rep)]
fwrite(jac, file.path(OUT, "subsampling_jaccard.csv"))
print(jac[, .(jaccard = mean(jaccard)), by = method][order(-jaccard)])

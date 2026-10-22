# Gene-list benchmark on GWAS Catalog disease gene lists.
#
# Query lists: genes mapped (MAPPED_GENE) to genome-wide significant
# associations (P <= 5e-8) for the diseases that have MalaCards KEGG relevance
# rankings in GSEABenchmarkeR.
#
# Gold standards
#   (1) MalaCards relevance (GSEABenchmarkeR evalRelevance, % of optimum)
#   (2) target-pathway rank: rank of the disease's own KEGG disease pathway
#   (3) temporal validation: methods are run on genes first reported up to the
#       disease's median discovery date ("early"); truth = KEGG pathways
#       over-represented (hypergeometric P < 0.01) among genes first reported
#       later and not in the early list ("late"). Scored by AUROC / AUPRC.
#
# Output: analysis/results/gwas/
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(GSEABenchmarkeR)
  library(org.Hs.eg.db)
  library(parallel)
})
OUT <- file.path(RES, "gwas")
dir.create(OUT, showWarnings = FALSE)

# ---------------------------------------------------------------- disease definitions
traits <- list(
  ALZ  = "^Alzheimer disease$",
  PARK = "^Parkinson disease$",
  CRC  = "^colorectal cancer$",
  BRCA = "^breast carcinoma$",
  PRAD = "^prostate carcinoma$",
  LES  = "^systemic lupus erythematosus$",
  DCM  = "^(idiopathic )?dilated cardiomyopathy$",
  PDCO = "^chronic obstructive pulmonary disease$",
  THCA = "^(differentiated )?thyroid (carcinoma|cancer)$",
  OV   = "^ovarian (serous )?carcinoma$",
  LAML = "^acute myeloid leukemia$",
  LGG  = "^glioma$",
  GBM  = "^glioblastoma$",
  PAAD = "^(exocrine )?pancreatic carcinoma$",
  KIRC = "^renal cell carcinoma$",
  UCEC = "^endometrial (carcinoma|cancer)$",
  LUAD = "^lung adenocarcinoma$",
  LUSC = "^squamous cell lung carcinoma$",
  SKCM = "^(cutaneous )?melanoma$",
  STAD = "^gastric (carcinoma|cancer)$",
  BLCA = "^urinary bladder (cancer|carcinoma)$",
  DMND = "^type 2 diabetes mellitus$",
  LIHC = "^hepatocellular carcinoma$",
  TGCT = "^testicular (germ cell tumor|cancer)$",
  CESC = "^cervical (carcinoma|cancer)$",
  DLBC = "^diffuse large B-cell lymphoma$",
  HNSC = "^head and neck (squamous cell carcinoma|cancer)$",
  UVM  = "^uveal melanoma$"
)
# KEGG disease pathways used as targets
targets <- c(ALZ = "hsa05010", PARK = "hsa05012", CRC = "hsa05210", BRCA = "hsa05224",
             PRAD = "hsa05215", LES = "hsa05322", DCM = "hsa05414", THCA = "hsa05216",
             LAML = "hsa05221", LGG = "hsa05214", GBM = "hsa05214", PAAD = "hsa05212",
             KIRC = "hsa05211", UCEC = "hsa05213", LUAD = "hsa05223", LUSC = "hsa05223",
             SKCM = "hsa05218", STAD = "hsa05226", BLCA = "hsa05219", DMND = "hsa04930",
             LIHC = "hsa05225", DLBC = NA, OV = NA, PDCO = NA, TGCT = NA, CESC = NA,
             HNSC = NA, UVM = NA)

nets <- readRDS(file.path(AN, "data", "nf_networks.rds"))
net <- nets$string_nokb400
gw <- readRDS(file.path(AN, "data", "gwas_sig.rds"))
split_genes <- function(x) {
  g <- unlist(strsplit(x, " - |, |; |,|;"))
  g <- trimws(g)
  g[nzchar(g)]
}
lists <- lapply(names(traits), function(code) {
  a <- gw[grepl(traits[[code]], MAPPED_TRAIT, ignore.case = TRUE)]
  if (!nrow(a)) return(NULL)
  long <- a[, .(gene = split_genes(MAPPED_GENE)), by = .(DATE, PUBMEDID, row = seq_len(nrow(a)))]
  long <- long[gene %in% net$genes]
  first <- long[, .(first = min(as.IDate(DATE))), by = gene]
  list(code = code, genes = first$gene, first = first)
})
names(lists) <- names(traits)
lists <- Filter(function(x) !is.null(x) && length(x$genes) >= 15, lists)
cat("diseases with >= 15 network genes:", length(lists), "\n")
print(vapply(lists, function(x) length(x$genes), 1L))

# ---------------------------------------------------------------- gene sets (KEGG, hsa ids, symbols)
kegg.gs <- readRDS(file.path(AN, "data", "kegg_gs_entrez.rds"))
eg2sym <- AnnotationDbi::mapIds(org.Hs.eg.db, unique(unlist(kegg.gs)), "SYMBOL", "ENTREZID")
gs <- lapply(kegg.gs, function(s) unique(unname(stats::na.omit(eg2sym[s]))))

# temporal split
temporal <- lapply(lists, function(x) {
  cut <- stats::median(x$first$first)
  early <- x$first[first <= cut, gene]
  late <- x$first[first > cut, gene]
  list(cut = cut, early = early, late = setdiff(late, early))
})
temporal <- Filter(function(t) length(t$early) >= 15 && length(t$late) >= 15, temporal)

jobs <- c(lapply(names(lists), function(c) list(id = c, kind = "full", q = lists[[c]]$genes)),
          lapply(names(temporal), function(c) list(id = c, kind = "early", q = temporal[[c]]$early)))

f_raw <- file.path(OUT, "gwas_raw.rds")
METHODS <- nf_methods_env(c("ora", "nea", "rwr", "guney", "v1", "nfr"))
rerun <- nzchar(Sys.getenv("NF_METHODS"))
if (!file.exists(f_raw) || rerun) {
  cl <- makeCluster(as.integer(Sys.getenv("NF_CORES", 14)))
  clusterExport(cl, c("gs", "net", "METHODS"))
  clusterEvalQ(cl, {
    source("analysis/scripts/00_setup.R")
    ctx <- make_ctx(gs, net)
    NULL
  })
  clusterSetRNGStream(cl, 99)
  res <- parLapplyLB(cl, jobs, function(j) {
    r <- run_genelist_methods(j$q, ctx, methods = METHODS, timing = TRUE)
    tm <- attr(r, "timing")
    r[, `:=`(disease = j$id, kind = j$kind, n_query = sum(net$genes %in% j$q))]
    list(r = r, time = tm)
  })
  stopCluster(cl)
  pv_new <- rbindlist(lapply(res, `[[`, "r"))
  if (rerun && file.exists(f_raw)) pv_new <- merge_results(readRDS(f_raw)$pvalues, pv_new)
  saveRDS(list(pvalues = pv_new, lists = lists,
               temporal = temporal, time = lapply(res, `[[`, "time")), f_raw)
}
pv <- readRDS(f_raw)$pvalues

# ---------------------------------------------------------------- (1) MalaCards relevance
mala <- readRDS(system.file("extdata/malacards/KEGG.rds", package = "GSEABenchmarkeR"))
full <- pv[kind == "full"]
ea.ranks <- lapply(split(full, full$method), function(dm) {
  lapply(split(dm, dm$disease), function(dd) {
    dd <- dd[!is.na(p)][order(p)]
    S4Vectors::DataFrame(GENE.SET = dd$set, PVAL = dd$p)
  })
})
d2p <- stats::setNames(names(lists), names(lists))
rel <- evalRelevance(ea.ranks, mala, d2p)
saveRDS(rel, file.path(OUT, "relevance.rds"))
print(round(sort(colMeans(rel), decreasing = TRUE), 2))

# ---------------------------------------------------------------- (2) target pathway rank
tg <- full[, .(target = targets[disease]), by = disease][!is.na(target)]
trank <- full[!is.na(p), .(set_id = sub("_.*", "", set), rank = rank(p) / .N,
                           p = p, q = p.adjust(p, "BH")), by = .(method, disease)]
trank <- merge(trank, unique(tg), by = "disease")[set_id == target]
fwrite(trank, file.path(OUT, "target_rank.csv"))
print(trank[, .(median_rank = median(rank), n_fdr = sum(q < 0.05)), by = method][order(median_rank)])

# ---------------------------------------------------------------- (3) temporal validation
auc <- function(score, truth) {
  r <- rank(score)
  n1 <- sum(truth)
  n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aupr <- function(score, truth) {
  o <- order(score, decreasing = TRUE)
  t <- truth[o]
  prec <- cumsum(t) / seq_along(t)
  sum(prec[t]) / sum(t)
}
universe <- unique(unlist(gs))
tv <- rbindlist(lapply(names(temporal), function(c) {
  late <- intersect(temporal[[c]]$late, universe)
  ora_late <- nf_ora(late, gs, universe = universe, min_size = 1L)
  truth_sets <- ora_late$set[ora_late$pvalue < 0.01]
  e <- pv[kind == "early" & disease == c & !is.na(p)]
  e[, .(auroc = auc(-log10(p), set %in% truth_sets),
        auprc = aupr(-log10(p), set %in% truth_sets),
        n_truth = sum(set %in% truth_sets), n_sets = .N), by = method][, disease := c]
}))
tv <- tv[n_truth >= 3]
fwrite(tv, file.path(OUT, "temporal_validation.csv"))
print(tv[, .(auroc = mean(auroc), auprc = mean(auprc), n = .N), by = method][order(-auroc)])

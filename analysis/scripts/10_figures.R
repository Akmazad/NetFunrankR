# Manuscript figures. Each block runs only if its results exist.
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("R/plotFactory.R")   # latest plotting code without reinstalling
suppressPackageStartupMessages({
})
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
ink <- "grey20"
fam_scale <- scale_fill_manual(values = FAMILY_COLS, name = NULL)
fam_col <- scale_colour_manual(values = FAMILY_COLS, name = NULL)
ord <- function(x, lv = METHOD_LEVELS) factor(x, levels = intersect(c(lv, setdiff(unique(x), lv)), unique(x)))

# ================================================================ Fig 2: calibration
f <- file.path(RES, "calibration_pvalues.rds")
if (file.exists(f)) {
  cal <- readRDS(f)[!is.na(p)]
  cal[, method := ord(method)]
  cal[, scheme := factor(scheme, c("uniform", "hub"),
                         c("Uniform random lists", "Hub-biased random lists"))]
  # (a) QQ, n = 100
  qq <- cal[size == 100, {
    pp <- sort(p)
    k <- unique(round(exp(seq(0, log(.N), length.out = 400))))
    .(exp = -log10(k / (.N + 1)), obs = -log10(pp[k]))
  }, by = .(method, scheme)]
  qq[, family := method_family(as.character(method))]
  ymax <- 8
  ends <- qq[, .SD[which.max(exp)], by = .(method, scheme)][
    method %in% c("ORA", "NEA", "NFR-combined", "Proximity-z", "NetFunrankR-v1")]
  ends[, lab := as.character(method)]
  ends[, y := pmin(obs, ymax)]
  # off-scale curves are labelled where they leave the panel
  exits <- qq[obs > ymax, .SD[which.min(exp)], by = .(method, scheme)]
  ends <- ends[!exits, on = .(method, scheme)]
  ends <- rbind(ends, exits[, .(method, scheme, exp, obs, family,
                                lab = paste0(method, " (off scale)"), y = ymax,
                                off = TRUE)],
                fill = TRUE)
  ends[is.na(off), off := FALSE]
  ends[off == TRUE, `:=`(exp = exp + 0.06, y = ymax - 0.8 * (seq_len(.N) - 1)), by = scheme]
  pa <- ggplot(qq, aes(exp, obs, group = method, colour = family)) +
    geom_abline(colour = "grey60", linetype = "dashed") +
    geom_line(linewidth = 0.6) +
    facet_wrap(~scheme) + fam_col +
    geom_text(data = ends[off == FALSE], aes(y = y, label = lab), hjust = 1.02, vjust = -0.4,
              size = 2.6, colour = ink, check_overlap = TRUE) +
    geom_text(data = ends[off == TRUE], aes(y = y, label = lab), hjust = 0, vjust = 1.2,
              size = 2.6, colour = ink) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = expression(Expected~-log[10]~p), y = expression(Observed~-log[10]~p),
         title = "a  Null p-value distributions (lists of 100 genes)") +
    theme_nf() + theme(legend.position = "none")
  # (b) FPR at 0.05
  s <- cal[, .(fpr = mean(p < 0.05)), by = .(method, scheme, size)]
  s[, family := method_family(as.character(method))]
  pb <- ggplot(s, aes(factor(size), fpr, fill = family)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.75) +
    facet_grid(scheme ~ method) + fam_scale +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Query size", y = "False-positive rate (p < 0.05)",
         title = "b  Type-I error per gene set") +
    theme_nf() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  # (c) any FDR hit per null list
  a <- cal[, .(any = any(p.adjust(p, "BH") < 0.05)), by = .(method, scheme, size, rep)][
    , .(frac = mean(any)), by = .(method, scheme, size)]
  a[, family := method_family(as.character(method))]
  pc <- ggplot(a, aes(method, frac, fill = family)) +
    geom_col(width = 0.75) + facet_grid(scheme ~ size, labeller = labeller(size = function(x) paste0("n = ", x))) +
    fam_scale + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "Null lists with >= 1 FDR < 0.05 pathway",
         title = "c  Spurious discoveries from random gene lists") +
    theme_nf() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
  save_fig(pa / (pb / pc + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
             plot_layout(heights = c(1, 2.5)), "Fig2_calibration", 10, 13)
  fwrite(s, file.path(RES, "calibration_fpr_table.csv"))
}

# ================================================================ Fig 3: power under incomplete annotation
f <- file.path(RES, "power_holdout.rds")
if (file.exists(f)) {
  pw <- readRDS(f)
  power_fig <- function(key, name) {
    s <- pw[method %in% key, .(power = mean(q < 0.05, na.rm = TRUE),
                               rank = median(rank, na.rm = TRUE)), by = .(method, s, h)]
    s[, method := factor(method, key)]
    s[, main := method == "NFR-combined"]
    facet <- facet_wrap(~s, labeller = labeller(s = function(x) paste0("signal = ", as.numeric(x) * 100, "% of query")))
    base <- function(y, ylab, title) {
      ggplot(s, aes(h, .data[[y]], colour = method, linewidth = main)) +
        geom_line() + geom_point(size = 1.6) + facet +
        scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.5), guide = "none") +
        scale_colour_manual(values = LINE_COLS, name = NULL) +
        scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
        labs(x = "Fraction of pathway genes missing from the annotation", y = ylab, title = title) +
        theme_nf()
    }
    p <- base("power", "Power (target FDR < 0.05)", "a  Detection of the target pathway") /
      base("rank", "Median rank percentile of target\n(lower = better)", "b  Ranking of the target pathway") +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
    save_fig(p, name, 10, 8)
    s
  }
  power_fig(c("ORA", "NEA", "RWR", "Proximity-z", "NFR-overlap", "NFR-combined"),
            "Fig3_power_incomplete_annotation")
  s <- power_fig(c("NFR-overlap", "NFR-linker", "NFR-distance", "NFR-combined", "ORA"),
                 "FigS_power_components")
  fwrite(dcast(s, method + s ~ h, value.var = "power"), file.path(RES, "power_table.csv"))
}

# ================================================================ Fig 4: GSEABenchmarkeR
f <- file.path(RES, "gseabench", "relevance.rds")
if (file.exists(f)) {
  rel <- readRDS(f)
  nsig <- readRDS(file.path(RES, "gseabench", "nr_sig_sets.rds"))
  rt <- readRDS(file.path(RES, "gseabench", "runtimes.rds"))
  rl <- melt(as.data.table(rel, keep.rownames = "dataset"), id.vars = "dataset",
             variable.name = "method", value.name = "relevance")
  ebmap <- c(ora = "ORA", camera = "CAMERA", gsva = "GSVA", globaltest = "GlobalTest", safe = "SAFE",
             gsea = "GSEA", padog = "PADOG", roast = "ROAST", samgs = "SAMGS")
  fix4 <- function(dt) {
    dt[, method := as.character(method)]
    dt <- dt[!method %in% c("NFR-total", "NFR-distance", "NFR-combined-dist")]
    dt[method %in% names(ebmap), method := ebmap[method]]
    dt
  }
  rl <- fix4(rl)
  lv <- rl[, .(m = median(relevance, na.rm = TRUE)), by = method][order(m), method]
  rl[, method := factor(method, lv)]
  rl[, family := method_family(as.character(method))]
  p4a <- ggplot(rl, aes(relevance, method, fill = family)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, colour = "grey30", linewidth = 0.3) +
    fam_scale + labs(x = "Phenotype relevance (% of optimal score)", y = NULL,
                     title = "a  MalaCards relevance, 42 GEO datasets") + theme_nf()
  ns <- melt(as.data.table(nsig, keep.rownames = "dataset"), id.vars = "dataset",
             variable.name = "method", value.name = "pct")
  ns <- fix4(ns)
  ns[, method := factor(as.character(method), lv)]
  ns[, family := method_family(as.character(method))]
  p4b <- ggplot(ns, aes(pct, method, fill = family)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, colour = "grey30", linewidth = 0.3) +
    fam_scale + labs(x = "% of gene sets significant (BH < 0.05)", y = NULL,
                     title = "b  Fraction of significant sets") + theme_nf()
  ebt <- melt(as.data.table(rt$eb, keep.rownames = "dataset"), id.vars = "dataset",
              variable.name = "method", value.name = "sec")
  nt <- rbindlist(lapply(rt$net, function(x) as.data.table(as.list(x))), fill = TRUE)
  nt <- melt(nt, measure.vars = names(nt), variable.name = "method", value.name = "sec")
  nt[method == "NFR", method := "NFR-combined"]
  tt <- rbind(ebt[, .(method = as.character(method), sec)], nt[, .(method = as.character(method), sec)])
  tt <- fix4(tt)
  tt[, family := method_family(method)]
  tt[, method := factor(method, unique(tt[, .(m = median(sec, na.rm = TRUE)), by = method][order(m), method]))]
  p4c <- ggplot(tt, aes(sec, method, fill = family)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, colour = "grey30", linewidth = 0.3) +
    scale_x_log10() + fam_scale +
    labs(x = "Runtime per dataset (s, log scale)", y = NULL, title = "c  Runtime") + theme_nf()
  save_fig((p4a | p4b | p4c) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
           "Fig4_GSEABenchmarkeR", 13, 5.5)
  fwrite(rl[, .(mean = mean(relevance), median = median(relevance)), by = method][order(-median)],
         file.path(RES, "gseabench", "relevance_summary.csv"))
}

# ================================================================ Fig 5: GWAS
f <- file.path(RES, "gwas", "relevance.rds")
if (file.exists(f)) {
  main_m <- c("ORA", "NEA", "RWR", "Proximity-z", "NetFunrankR-v1", "NFR-overlap", "NFR-combined")
  box_panel <- function(dt, y, ylab, title, methods, hline = NULL, pct = FALSE) {
    dt <- dt[method %in% methods]
    dt[, method := factor(method, methods)]
    dt[, family := method_family(as.character(method))]
    p <- ggplot(dt, aes(method, .data[[y]], fill = family))
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey40")
    p <- p + geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey30", linewidth = 0.3) +
      geom_point(position = position_jitter(width = 0.12, seed = 1), size = 0.7, colour = "grey25") +
      fam_scale + labs(x = NULL, y = ylab, title = title) +
      theme_nf() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (pct) p <- p + scale_y_continuous(labels = scales::percent)
    p
  }
  rel <- readRDS(f)
  rl <- melt(as.data.table(rel, keep.rownames = "disease"), id.vars = "disease",
             variable.name = "method", value.name = "relevance")
  rl[, method := as.character(method)]
  p5a <- box_panel(rl, "relevance", "Relevance (% of optimum)",
                   "a  MalaCards relevance (25 diseases)", main_m)
  # paired gains over ORA
  w <- dcast(rl, disease ~ method, value.var = "relevance")
  gm <- c("NEA", "RWR", "Proximity-z", "NFR-overlap", "NFR-combined")
  gain <- rbindlist(lapply(gm, function(m) data.table(method = m, disease = w$disease,
                                                       gain = w[[m]] - w$ORA)))
  pv <- gain[, .(p = wilcox.test(gain)$p.value, wins = sum(gain > 0), n = .N,
                 med = median(gain)), by = method]
  pv[, lab := sprintf("%d/%d, p=%s", wins, n, formatC(p, format = "g", digits = 2))]
  gain[, method := factor(method, gm)]
  pv[, method := factor(method, gm)]
  gain[, family := method_family(as.character(method))]
  p5b <- ggplot(gain, aes(gain, method, colour = family)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_point(position = position_jitter(height = 0.15, seed = 2), size = 1.2, alpha = 0.8) +
    geom_point(data = pv, aes(x = med, y = method), inherit.aes = FALSE, shape = 124, size = 6) +
    geom_text(data = pv, aes(x = max(gain$gain) + 1, y = method, label = lab), inherit.aes = FALSE,
              hjust = 0, size = 2.7, colour = ink) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.45))) + fam_col +
    labs(x = "Relevance gain over ORA (percentage points, per disease)", y = NULL,
         title = "b  Paired gain over ORA (| = median)") + theme_nf() + guides(colour = "none")
  # significant calls
  sc <- fread(file.path(RES, "gwas", "significant_call_precision.csv"))[method %in% main_m]
  sc <- melt(sc[, .(method, Relevant = n_sig_relevant, Other = n_sig - n_sig_relevant, precision)],
             id.vars = c("method", "precision"), variable.name = "kind", value.name = "n")
  sc[, method := factor(method, rev(main_m))]
  sc[, kind := factor(kind, c("Other", "Relevant"))]
  p5c <- ggplot(sc, aes(n, method, fill = kind)) +
    geom_col(width = 0.65) +
    geom_text(data = sc[, .(n = sum(n), lab = sprintf("%d/%d (%.0f%%)", n[kind == "Relevant"], sum(n), 100 * precision[1])), by = method],
              aes(x = n, y = method, label = lab), inherit.aes = FALSE, hjust = -0.05, size = 2.6, colour = ink) +
    scale_fill_manual(values = c(Relevant = "#eb6834", Other = "grey78"), name = NULL,
                      breaks = c("Relevant", "Other")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = "Pathways called at FDR < 0.05 (summed over 25 diseases)", y = NULL,
         title = "c  Significant calls: MalaCards-relevant vs other") + theme_nf()
  tr <- fread(file.path(RES, "gwas", "target_rank.csv"))
  p5d <- box_panel(tr, "rank", "Rank percentile of KEGG disease\npathway (lower = better)",
                   "d  Target disease pathway (21 diseases)", main_m, pct = TRUE)
  tv <- fread(file.path(RES, "gwas", "temporal_validation.csv"))
  p5e <- box_panel(tv, "auroc", "AUROC", "e  Temporal validation (10 diseases)", main_m, hline = 0.5)
  save_fig(((p5a | p5b) / (p5c | p5d | p5e)) + plot_layout(guides = "collect") &
             theme(legend.position = "bottom"), "Fig5_GWAS_benchmark", 14, 10)
  # supplement: all NetFunrankR components
  comp <- c("ORA", "NFR-overlap", "NFR-linker", "NFR-distance", "NFR-total", "NFR-combined-dist", "NFR-combined")
  save_fig((box_panel(rl, "relevance", "Relevance (% of optimum)", "a  MalaCards relevance", comp) |
              box_panel(tr, "rank", "Target pathway rank percentile", "b  Target pathway", comp, pct = TRUE) |
              box_panel(tv, "auroc", "AUROC", "c  Temporal validation", comp, hline = 0.5)) +
             plot_layout(guides = "collect") & theme(legend.position = "bottom"),
           "FigS_GWAS_components", 14, 5)
}

# ================================================================ Fig 6: robustness
f <- file.path(RES, "robustness", "robustness_relevance.csv")
if (file.exists(f)) {
  rr <- fread(f)
  s <- rr[, .(rel = mean(rel)), by = .(variant, rep, method)][
    , .(rel = mean(rel), lo = min(rel), hi = max(rel)), by = .(variant, method)]
  vlev <- c("STRING noKB>=400", "STRING noKB>=700", "STRING all>=700 (circular)",
            "bin size 50", "bin size 250", "dropout 10%", "dropout 25%", "dropout 50%", "rewired")
  s <- s[method %in% c("ORA", "Proximity-z", "NFR-overlap", "NFR-linker", "NFR-distance", "NFR-combined")]
  s[, variant := factor(variant, rev(vlev))]
  s[, method := ord(method)]
  p6a <- ggplot(s, aes(rel, variant, colour = method)) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
    scale_colour_manual(values = LINE_COLS, name = NULL) +
    labs(x = "Mean MalaCards relevance of GWAS lists (%)", y = NULL,
         title = "a  Interactome variants, perturbations and rewired control") + theme_nf()
  jac <- fread(file.path(RES, "robustness", "subsampling_jaccard.csv"))
  jac[, method := ord(method)]
  jac[, family := method_family(as.character(method))]
  p6b <- ggplot(jac, aes(method, jaccard, fill = family)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6, colour = "grey30", linewidth = 0.3) +
    fam_scale + labs(x = NULL, y = "Jaccard (top-20, half vs full list)",
                     title = "b  Stability to query subsampling") +
    theme_nf() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p6a | p6b, "Fig6_robustness", 12, 5)
}

# ================================================================ Fig 7: COPD case study
f <- file.path(RES, "copd", "copd_kegg_netfunrankr_v2.csv")
if (file.exists(f)) {
  net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
  res <- fread(f)
  res[, set := sub("^hsa[0-9]+_", "", set)][, set := gsub("_", " ", set)]
  class(res) <- c("nf_result", class(res))
  p7a <- plot_evidence_map(res, n_label = 14) + ggtitle("a  Evidence map, COPD genes vs KEGG")
  kegg.gs <- readRDS(file.path(AN, "data", "kegg_gs_entrez.rds"))
  eg2sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(unlist(kegg.gs)),
                                  "SYMBOL", "ENTREZID")
  raw <- fread(f)
  netonly <- raw[evidence %in% c("network only", "overlap + network")][order(p_beyond)]
  if (nrow(netonly)) {
    s1 <- netonly$set[1]
    Y <- unique(stats::na.omit(eg2sym[kegg.gs[[s1]]]))
    p7b <- plot_linker_network(COPD_proteinCodingGenes, Y, net, max_members = 35) +
      ggtitle(paste0("b  ", gsub("_", " ", sub("^hsa[0-9]+_", "", s1))))
  } else p7b <- plot_spacer()
  sm <- fread(file.path(RES, "copd", "copd_small_lists.csv"))
  s <- sm[, .(rel = mean(relevance), se = sd(relevance) / sqrt(.N)), by = .(method, size)]
  s <- s[method %in% c("ORA", "NEA", "RWR", "Proximity-z", "NFR-overlap", "NFR-combined")]
  s[, method := factor(method, names(LINE_COLS))]
  p7c <- ggplot(s, aes(size, rel, colour = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = rel - se, ymax = rel + se), width = 1.2) +
    scale_colour_manual(values = LINE_COLS, name = NULL) +
    labs(x = "Number of COPD genes in query", y = "MalaCards COPD relevance (%)",
         title = "c  Small query lists (50 random subsets each)") + theme_nf()
  save_fig((p7a | p7b) / p7c + plot_layout(heights = c(1.3, 1)), "Fig7_COPD_case", 13, 10)
}

# ================================================================ Fig S: runtime
f <- file.path(RES, "runtime.csv")
if (file.exists(f)) {
  rt0 <- fread(f)
  v1 <- rt0[grepl("original", method)]
  rt <- rt0[!is.na(size) & !is.na(seconds)]
  # m_nfr() runs the test twice (linker and distance combinations)
  rt[method == "NFR (all four tests)", `:=`(seconds = seconds / 2, method = "NetFunrankR v2")]
  rt[, method := sub(" [(].*", "", method)]
  rt[method == "NetFunrankR-v1 statistic", method := "NetFunrankR v1 statistic, 1000 perm"]
  rt[method == "Proximity-z", method := "Proximity-z, 1000 perm"]
  rt[method == "RWR", method := "RWR, 200 perm"]
  rt[, seconds := pmax(seconds, 0.01)]          # timer resolution
  cols <- c(`NetFunrankR v2` = "#eb6834", ORA = "#2a78d6", `Proximity-z, 1000 perm` = "#1baf7a",
            `RWR, 200 perm` = "#eda100", `NetFunrankR v1 statistic, 1000 perm` = "#e87ba4",
            NEA = "#6250d6")
  p8 <- ggplot(rt, aes(size, seconds, colour = method)) +
    geom_line(aes(linewidth = method == "NetFunrankR v2")) + geom_point(size = 1.6) +
    scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.4), guide = "none") +
    facet_wrap(~library) + scale_x_log10() +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                  labels = c("<0.01", "0.1", "1", "10", "100", "1,000", "10,000")) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = "Query genes (log scale)", y = "Seconds for the whole library (log scale)") +
    theme_nf()
  if (nrow(v1)) {
    lab <- data.table(library = "KEGG", size = v1$size, seconds = v1$seconds,
                      txt = sprintf("NetFunrankR v1 original code\n(100 perm): %.1f h", v1$seconds / 3600))
    p8 <- p8 + geom_point(data = lab, aes(size, seconds), inherit.aes = FALSE, shape = 4, size = 3) +
      geom_text(data = lab, aes(size, seconds, label = txt), inherit.aes = FALSE,
                hjust = 1.08, size = 2.8, colour = ink)
  }
  save_fig(p8, "FigS_runtime", 10, 4.8)
}

# Accuracy of the exact (convolution) and saddlepoint null distributions.
#  (a) exact vs Monte-Carlo permutation (2,000 degree-binned permutations),
#      total proximity test, COPD genes vs KEGG
#  (b) exact vs saddlepoint, random hub-biased lists of 100 and 400 genes
# Output: analysis/results/null_accuracy.csv, figures/FigS_null_accuracy.*
source("analysis/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
gs <- kegg_symbols()
D <- geneSetDistances(net, gs)

set.seed(21)
ex <- nfTest(COPD_proteinCodingGenes, gs, net, D = D)
mc <- nfPerm(COPD_proteinCodingGenes, gs, net, D = D, n_perm = 2000L)
a <- merge(as.data.table(ex)[, .(set, exact = p_proximity)],
           as.data.table(mc)[, .(set, mc = p_perm)], by = "set")
cat("exact vs MC Spearman (log p):", cor(log10(a$exact), log10(a$mc), method = "spearman"), "\n")

qs <- list(`random, 100 genes` = sample(nodeNames(net), 100, prob = nodeDegree(net)),
           `random, 400 genes` = sample(nodeNames(net), 400, prob = nodeDegree(net)),
           `COPD, 261 genes` = COPD_proteinCodingGenes)
b <- rbindlist(lapply(names(qs), function(n) {
  q <- qs[[n]]
  e <- as.data.table(nfTest(q, gs, net, D = D, exact_max = 1e6))
  s <- as.data.table(nfTest(q, gs, net, D = D, exact_max = 0))
  cols <- c("p_overlap", "p_linker", "p_beyond", "p_proximity")
  m <- merge(e[, c("set", cols), with = FALSE], s[, c("set", cols), with = FALSE],
             by = "set", suffixes = c(".exact", ".saddle"))
  rbindlist(lapply(cols, function(cl)
    data.table(n = n, test = sub("p_", "", cl), exact = m[[paste0(cl, ".exact")]],
               saddle = m[[paste0(cl, ".saddle")]])))
}))
b[, dlog := abs(log10(exact) - log10(saddle))]
print(b[, .(max_dlog10 = max(dlog), max_dlog10_p_lt_1e3 = max(c(0, dlog[exact < 1e-3]))), by = .(n, test)])
fwrite(rbind(a[, .(panel = "a", test = "proximity", exact, other = mc)],
             b[, .(panel = paste0("b: ", n), test, exact, other = saddle)]),
       file.path(RES, "null_accuracy.csv"))

floor_mc <- -log10(1 / 2001)
pa <- ggplot(a, aes(-log10(exact), -log10(mc))) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = floor_mc, linetype = "dotted", colour = "grey50") +
  geom_point(size = 1.4, alpha = 0.7, colour = "#2a78d6") +
  annotate("text", x = 0.2, y = floor_mc, label = "permutation floor 1/(B+1)",
           vjust = -0.5, hjust = 0, size = 3, colour = "grey30") +
  labs(x = expression(-log[10]~p~"(exact)"), y = expression(-log[10]~p~"(2,000 permutations)"),
       title = "a  Exact vs permutation") + theme_nf()
b[, test := factor(test, c("overlap", "linker", "beyond", "proximity"),
                   c("overlap", "linker", "distance", "total proximity"))]
pb <- ggplot(b, aes(-log10(exact), -log10(saddle), colour = test)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~factor(n, names(qs)), scales = "free") +
  scale_colour_manual(values = c(overlap = "#2a78d6", linker = "#eb6834",
                                 distance = "#1baf7a", `total proximity` = "#eda100"), name = NULL) +
  labs(x = expression(-log[10]~p~"(exact)"), y = expression(-log[10]~p~"(saddlepoint)"),
       title = "b  Exact vs saddlepoint") + theme_nf()
p <- pa | pb
p <- p + plot_layout(widths = c(1, 3))
ggsave(file.path(FIG, "FigS_null_accuracy.pdf"), p, width = 14, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIG, "FigS_null_accuracy.png"), p, width = 14, height = 4.5, dpi = 200, bg = "white")

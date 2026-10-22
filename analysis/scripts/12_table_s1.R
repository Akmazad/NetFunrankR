# Supplementary Table S1: GWAS gene lists
suppressPackageStartupMessages(library(data.table))
raw <- readRDS("analysis/results/gwas/gwas_raw.rds")
src <- readLines("analysis/scripts/06_gwas_benchmark.R")
i1 <- grep("^traits <- list", src)
i2 <- i1 + which(src[(i1 + 1):length(src)] == ")")[1]
traits <- eval(parse(text = paste(src[i1:i2], collapse = "\n")))
j1 <- grep("^targets <- c", src)
j2 <- j1 + which(grepl("\\)$", src[j1:length(src)]))[1] - 1
targets <- eval(parse(text = paste(src[j1:j2], collapse = "\n")))
tab <- rbindlist(lapply(names(raw$lists), function(c) {
  t <- raw$temporal[[c]]
  data.table(Code = c,
             Trait = gsub("[$^]", "", traits[[c]]),
             Genes = length(raw$lists[[c]]$genes),
             Target = ifelse(is.na(targets[c]), "--", targets[c]),
             Split = if (is.null(t)) "--" else format(t$cut),
             Early = if (is.null(t)) "--" else as.character(length(t$early)),
             Late = if (is.null(t)) "--" else as.character(length(t$late)))
}))
esc <- function(x) gsub("([&%_#])", "\\\\\\1", x)
lines <- c("\\begin{table}[h]\\centering\\small",
           paste0("\\caption{\\textbf{GWAS gene lists.} Genes mapped to genome-wide significant ",
                  "associations ($p\\le5\\times10^{-8}$) of the GWAS Catalog trait (EFO label pattern), ",
                  "restricted to the interactome; KEGG target disease pathway; temporal split date and ",
                  "numbers of early and late genes (temporal validation used diseases with $\\ge$15 genes ",
                  "in both parts and $\\ge$3 late-enriched pathways).}"),
           "\\label{tab:s1}",
           "\\begin{tabular}{llrllrr}\\toprule",
           "Code & GWAS trait & Genes & KEGG target & Split date & Early & Late\\\\\\midrule",
           tab[, sprintf("%s & %s & %d & %s & %s & %s & %s\\\\", Code, esc(Trait), Genes,
                         Target, Split, Early, Late)],
           "\\bottomrule\\end{tabular}\\end{table}")
writeLines(lines, "analysis/manuscript/tableS1.tex")
print(tab)

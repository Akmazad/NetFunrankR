
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NetFunrankR: an R package for robust functional enrichment test using Network medicine approach and beyond

## Introduction

With the given biomarker sets, this package primarily facilitates the
following analysis: 1. **Topology-aware enrichment**: Network-proximity
based functional enrichment test utilizing protein-protein interaction
network 2. **Over-representation enrichment** Typical
hypergeometric-distribution based over-representation test

Results include two spreadsheets, one for ORA-based test, and another
for TW-based test. Visualization includes, a venn diagram, double-bar
plots and two lollipop plots. This document shows the basic usage of
this package:

## Installation

``` r
BiocManager::install("EnsDb.Hsapiens.v75")
library(devtools)
install_github("Akmazad/NetFunrankR")
```

### Load Package

``` r
library(dplyr)
library(data.table)
library(ggplot2)
library(NetFunrankR)
```

### Load example disease genes

The package comes with an example vector of disease genes (COPD) named
example_enrich.df. Let’s observe the head of this data..

``` r
print(head(COPD_proteinCodingGenes))
```

### Perform analysis

With the given disease genes, a wrapper function (`NetFunrankR_wrapper`)
can be invoked to run the whole analysis, and the results (reports,
plots - both pdf and ggplot2 objects) are saved in a given directory
(default: `./NetFunrankR_Results/`).

``` r
NetFunrankR_wrapper(
  disease_genes = COPD_proteinCodingGenes,
  PathwayRanker_Type = "Both",
  nPermute = 10,
  outdir = "./NetFunrankR_Results/"
)
```

### Check results

#### Number of significantly enriched terms/pathways

This venn diagram shows counts of significantly enriched terms/pathways
by the two methods

``` r
load("./NetFunrankR_results/CommonFunctions_ORA_TW.vennplt.obj.RData")
common.terms.vennplt.obj
```

#### Significantly enriched terms/pathways

This plot shows the terms/pathways that are significantly enriched
(`pvalue < 0.05`) and commonly found between Network-Proximity based and
Over-representation based methods

``` r
load("./NetFunrankR_results/CommonFunctions_ORA_TW.barplt.obj.RData")
common.terms.barplt.obj
```

#### Exclusive terms/pathways found by ORA

This plot shows the terms/pathways that are exclusively enriched
(`pvalue < 0.05`) via ORA method

``` r
load("./NetFunrankR_results/Only_ORA_terms.plt.obj.RData")
only_ORA_terms.plt.obj
```

#### Exclusive terms/pathways found by TW

This plot shows the terms/pathways that are exclusively enriched
(`pvalue < 0.05`) via TW method

``` r
load("./NetFunrankR_results/Only_TW_terms.plt.obj.RData")
only_TW_terms.plt.obj
```

### Further configuration of results

The enrichment results from both the methods can be further controlled
based on filtering `p.adjust < 0.05` instead of the nominal `pvalue`.

``` r
    rankedPathways_ORA <- fread("NetFunrankR_results/KEGG_PE_ORA_result.csv")
    rankedPathways_TW <- fread("NetFunrankR_results/KEGG_PE_TopoAw_result.csv")
    rankedPathways_ORA_filt <- rankedPathways_ORA %>% dplyr::filter(p.adjust < 0.05)
    rankedPathways_TW_filt <- rankedPathways_TW  %>% dplyr::filter(p.adjust < 0.05)
```

Then the results can be visualized accordingly.

``` r
library(ggvenn)
library(tidyr)
library(forcats)
library(ggplot2)
library(cowplot)

# venn diagram
p1 <- plot_venn(
        vec1 = rankedPathways_ORA_filt %>% dplyr::pull(Description),
        vec2 = rankedPathways_TW_filt %>% dplyr::pull(Description))
# double-bar plot
p2 <- dplyr::inner_join(
        rankedPathways_ORA_filt,
        rankedPathways_TW_filt,
        by = "Description") %>%
  plot_double_bar()
# exclusive ORA 
p3 <- dplyr::anti_join(
        rankedPathways_ORA_filt,
        rankedPathways_TW_filt,
        by = c("Description" = "Description")) %>%
  plot_lollipop(color_palette = c("#616530FF", "#dce58f"))
# exclusive TW
p4 <- dplyr::anti_join(
        rankedPathways_TW_filt,
        rankedPathways_ORA_filt,
        by = c("Description" = "Description")) %>%
  plot_lollipop(color_palette = c("#CC8214FF", "#f6e8b1"))
cowplot::plot_grid(p1,p2,p3,p4,
                   # rel_heights = c(1,1),
                   # rel_widths = c(1,1), 
                   labels = c("a", "b", "c", "d"), 
                   ncol = 2,
                   align = "l"
                   )
```

### Performance analysis on smaller input

With the given disease genes, a wrapper function (`NetFunrankR_wrapper`)
can be invoked to run the whole analysis, and the results (reports,
plots - both pdf and ggplot2 objects) are saved in a given directory
(default: `./NetFunrankR_Results_small/`).

``` r
NetFunrankR_wrapper(
  disease_genes = COPD_proteinCodingGenes %>% head(6),
  PathwayRanker_Type = "Both",
  nPermute = 10,
  outdir = "./NetFunrankR_Results_small/"
)
```

# Conclusion

In this page, we demonstrated the basic functionality of the NetFunrankR
package. We loaded example disease genes, performed both
Over-representation (ORA) and topology-aware (TW) network-proximity
methods for functional term/pathway enrichment, compare, report and
visualize results. For more detailed analysis, refer to the package
documentation.

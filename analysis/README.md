# Reproducing the NetFunrankR manuscript

Scripts that reproduce every number, figure and table of the manuscript.
They are not part of the R package (excluded from the build).

## Run

From the repository root:

```bash
RSCRIPT=Rscript bash analysis/run_all.sh                 # everything
RSCRIPT=Rscript bash analysis/run_all.sh 06_gwas_benchmark 10_figures   # selected steps
```

`RSCRIPT` is the path to `Rscript` if it is not on the `PATH`. The scripts
install the package from this repository into a private library
(`tools::R_user_dir("NetFunrankR", "cache")`) so that results always come
from the code in the working tree. Requirements: R >= 4.4 with Bioconductor
(GSEABenchmarkeR, EnrichmentBrowser, org.Hs.eg.db), the CRAN packages
data.table, patchwork, ggplot2 and magick, and a C++ toolchain. All input
data are downloaded automatically. A complete run takes about 4-5 hours on a
20-core laptop.

## Steps

| Script | Output |
|---|---|
| `01_build_networks.R` | STRING download; noKB >= 400 (primary), noKB >= 700, full >= 700 networks |
| `02_prepare_inputs.R` | network objects, KEGG gene sets, GEO2KEGG compendium with limma DE, GWAS Catalog hits |
| `03_calibration.R` | type-I error under uniform and hub-biased random lists (Figure 2) |
| `04_power_annotation_holdout.R` | power under incomplete annotation (Figure 3) |
| `05_gseabenchmarker.R` | GSEABenchmarkeR, 42 expression datasets (Figure 4) |
| `06_gwas_benchmark.R` | GWAS gene lists: MalaCards relevance, target pathway, temporal validation (Figure 5) |
| `07_robustness.R` | interactome variants, edge dropout, rewiring, bin size, subsampling (Figure 6) |
| `08_runtime.R` | runtime, including the original permutation implementation (Figure S2) |
| `09_case_copd.R` | COPD case study (Figure 7) |
| `11_null_accuracy.R` | exact vs permutation vs saddlepoint (Figure S1) |
| `12_table_s1.R` | Supplementary Table S1 |
| `10_figures.R` | all figures |

## Settings (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `NF_CORES` | all cores but one | worker processes |
| `NF_LIB` | user cache directory | private package library |
| `NF_REINSTALL` | `FALSE` | reinstall the package before running |
| `NF_METHODS` | all | comma-separated subset of methods to (re)run, merged into existing results: `ora,nea,rwr,closest,symmetric,nfr` |
| `NF_STRING_VERSION`, `NF_TAXON` | `12.0`, `9606` | STRING release and organism |
| `NF_GWAS_URL`, `NF_GWAS_P` | latest release, `5e-8` | GWAS Catalog archive and significance threshold |
| `NF_NREP`, `NF_NTARGET` | 200, 100 | replicates (calibration), target pathways (power) |
| `NF_ORIGINAL_REF` | `169438c` | git revision with the original implementation (runtime baseline) |

## Method names

| Label | Method |
|---|---|
| ORA | hypergeometric over-representation analysis |
| NEA | network enrichment analysis (edge count, configuration-model z-score) |
| RWR | random walk with restart, degree-binned permutation null |
| Proximity-closest | closest-distance proximity d_c (Guney et al. 2016), permutation z-score |
| Proximity-symmetric | symmetric module distance d_XY (Menche et al. 2015), permutation z-score |
| **NetFunrankR** | combined test (degree-aware overlap + linker, Fisher) |
| NetFunrankR-overlap / -linker / -distance / -total | individual NetFunrankR tests |
| NetFunrankR-distance-combined | combination with the distance instead of the linker test |

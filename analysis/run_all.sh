#!/usr/bin/env bash
# Reproduce all benchmarks, figures and tables (from the package root):
#   bash analysis/run_all.sh                # everything (~4-5 h on a 20-core laptop)
#   bash analysis/run_all.sh 06_gwas_benchmark 10_figures   # selected steps
# Requirements: R >= 4.4 with Bioconductor (GSEABenchmarkeR, EnrichmentBrowser,
# org.Hs.eg.db), Rtools on Windows. Raw downloads go to analysis/data/raw:
#   STRING v12: 9606.protein.links.detailed.v12.0.txt.gz, 9606.protein.info.v12.0.txt.gz
#   GWAS Catalog: gwas-catalog-associations_ontology-annotated-full.zip (unzipped)
# Set NF_CORES to the number of worker processes (default: all cores but one).
RS="${RSCRIPT:-Rscript}"   # path to Rscript if it is not on the PATH
steps=${@:-01_build_networks 02_prepare_inputs 03_calibration 04_power_annotation_holdout 06_gwas_benchmark 07_robustness 08_runtime 09_case_copd 05_gseabenchmarker 11_null_accuracy 12_table_s1 10_figures}
mkdir -p analysis/results
NF_REINSTALL=TRUE "$RS" -e 'source("analysis/scripts/00_setup.R")' > analysis/results/log_install.txt 2>&1
for s in $steps; do
  echo "=== $s $(date)"
  "$RS" "analysis/scripts/$s.R" > "analysis/results/log_$s.txt" 2>&1
  echo "    exit $? $(date)"
done

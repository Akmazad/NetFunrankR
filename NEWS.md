# NetFunrankR 0.99.0

* Submission to Bioconductor.
* New decomposed network-proximity test (`netEnrich()`): degree-aware
  overlap test, beyond-overlap linker and distance tests conditional on the
  overlap, and their Fisher combination, all with exact null distributions
  (discrete convolution; saddlepoint approximation for long lists).
* S4 classes `NFNetwork` and `NFResult` with accessors and `show()` methods.
* Bundled STRING v12 interactome with a combined score recomputed without
  the curated-database and text-mining channels (`STRING_v12_PPI`,
  `stringNetwork()`).
* Gene sets can be supplied as lists, tables, GMT files or
  `GSEABase::GeneSetCollection` objects.
* `netEnrichBatch()` tests many gene lists in parallel with BiocParallel.
* Plots: `plotEvidenceMap()`, `plotTopSets()`, `plotLinkerNetwork()`,
  `plotNullDistribution()`.
* Interactive Shiny dashboard: `runNetFunrankRApp()`, configurable with any
  gene-set libraries, interactomes and example list (`dashboardDefaults()`),
  asynchronous computation with `mirai`, and deployment recipes (Docker,
  Hugging Face Spaces, Posit Connect Cloud, shinyapps.io) in `deploy/`.
* `geneSetDistances()` processes gene sets in blocks (`chunkSize`) to limit
  peak memory.
* Fixed in the previous (unreleased) version: off-by-one in the
  hypergeometric test, degree-matched sampling of unique-degree hubs, and
  circular use of pathway-derived interaction evidence.

# Runtime (single core).
#  (a) the original permutation-based implementation of symmetric proximity
#      (per-pathway network rebuild, exact-degree permutation), taken from the
#      git history, on 3 KEGG pathways and extrapolated to the library
#  (b) all methods across query sizes for KEGG and Reactome
# Output: analysis/results/runtime.csv
source("analysis/scripts/00_setup.R")

# git revision that contains the original implementation (NF_ORIGINAL_REF)
ORIGINAL_REF <- Sys.getenv("NF_ORIGINAL_REF", "169438c")
ORIGINAL_LABEL <- "Original implementation (symmetric proximity, 100 perm)"

# ---------------------------------------------------------------- (a) original
orig <- new.env()
tmp <- tempfile(fileext = ".rda")
system2("git", c("show", paste0(ORIGINAL_REF, ":data/protein.links.rda")), stdout = tmp)
load(tmp, envir = orig)
src <- tempfile(fileext = ".R")
system2("git", c("show", paste0(ORIGINAL_REF, ":R/netProximityMeasure.R")), stdout = src)
suppressPackageStartupMessages({
  library(dplyr)
})
sys.source(src, envir = orig)
orig$data <- function(...) invisible(NULL)        # protein.links already in env
gs_kegg <- kegg_symbols()
q <- COPD_proteinCodingGenes
t_orig <- vapply(names(gs_kegg)[1:3], function(s) {
  system.time(orig$quantify_network_proximity_score(900, q, gs_kegg[[s]], nPermute = 100))[3]
}, numeric(1))
orig_time <- data.table(method = ORIGINAL_LABEL, library = "KEGG",
                      size = length(q), n_perm = 100,
                      seconds = mean(t_orig) * length(gs_kegg))
cat("original implementation, per pathway (100 perm):", mean(t_orig), "s; extrapolated:", orig_time$seconds, "s\n")

# ---------------------------------------------------------------- (b) all methods
net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
libs <- list(KEGG = gs_kegg, Reactome = loadPathways("Reactome"))
sizes <- c(10L, 50L, 100L, 250L, 500L, 1000L, 2000L)
set.seed(11)
rt <- list()
for (lib in names(libs)) {
  t_ctx <- system.time(ctx <- make_ctx(libs[[lib]], net))[3]
  rt[[length(rt) + 1]] <- data.table(method = "set distances (once per library)",
                                     library = lib, size = NA, n_perm = NA, seconds = t_ctx)
  for (n in sizes) {
    qq <- sample(nodeNames(net), n, prob = nodeDegree(net))
    tm <- c(
      `NetFunrankR (all tests)` = system.time(m_nfr(qq, ctx))[3],
      `ORA` = system.time(m_ora(qq, ctx))[3],
      `NEA` = system.time(m_nea(qq, ctx))[3],
      `RWR (200 perm)` = system.time(m_rwr(qq, ctx))[3],
      `Proximity-closest (1000 perm)` = system.time(m_prox_closest(qq, ctx))[3],
      `Proximity-symmetric (1000 perm)` = if (n <= 500 || lib == "KEGG")
        system.time(m_prox_symmetric(qq, ctx, B = 1000L))[3] else NA_real_)
    rt[[length(rt) + 1]] <- data.table(method = sub("\\.elapsed$", "", names(tm)),
                                       library = lib, size = n,
                                       n_perm = NA, seconds = unname(tm))
    print(rt[[length(rt)]])
  }
}
rt <- rbindlist(c(list(orig_time), rt), fill = TRUE)
fwrite(rt, file.path(RES, "runtime.csv"))

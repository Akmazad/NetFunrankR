# Runtime (single core).
#  (a) original NetFunrankR v1 implementation (git HEAD~: per-pathway network
#      rebuild, exact-degree permutation) on 3 KEGG pathways, extrapolated
#  (b) v2 methods across query sizes for KEGG and Reactome
# Output: analysis/results/runtime.csv
source("analysis/scripts/00_setup.R")

# ---------------------------------------------------------------- (a) v1
v1 <- new.env()
tmp <- tempfile(fileext = ".rda")
system2("git", c("show", "169438c:data/protein.links.rda"), stdout = tmp)
load(tmp, envir = v1)
src <- tempfile(fileext = ".R")
system2("git", c("show", "169438c:R/netProximityMeasure.R"), stdout = src)
suppressPackageStartupMessages({
  library(dplyr)
})
sys.source(src, envir = v1)
v1$data <- function(...) invisible(NULL)        # protein.links already in env
gs_kegg <- kegg_symbols()
q <- COPD_proteinCodingGenes
t_v1 <- vapply(names(gs_kegg)[1:3], function(s) {
  system.time(v1$quantify_network_proximity_score(900, q, gs_kegg[[s]], nPermute = 100))[3]
}, numeric(1))
v1_time <- data.table(method = "NetFunrankR v1 (original code)", library = "KEGG",
                      size = length(q), n_perm = 100,
                      seconds = mean(t_v1) * length(gs_kegg))
cat("v1 per pathway (100 perm):", mean(t_v1), "s; extrapolated:", v1_time$seconds, "s\n")

# ---------------------------------------------------------------- (b) v2
net <- readRDS(file.path(AN, "data", "nf_networks.rds"))$string_nokb400
libs <- list(KEGG = gs_kegg, Reactome = load_pathways("Reactome"))
sizes <- c(10L, 50L, 100L, 250L, 500L, 1000L, 2000L)
set.seed(11)
rt <- list()
for (lib in names(libs)) {
  t_ctx <- system.time(ctx <- make_ctx(libs[[lib]], net))[3]
  rt[[length(rt) + 1]] <- data.table(method = "set distances (once per library)",
                                     library = lib, size = NA, n_perm = NA, seconds = t_ctx)
  for (n in sizes) {
    qq <- sample(net$genes, n, prob = net$degree)
    tm <- c(
      `NFR (all four tests)` = system.time(m_nfr(qq, ctx))[3],
      `ORA` = system.time(m_ora(qq, ctx))[3],
      `NEA` = system.time(m_nea(qq, ctx))[3],
      `RWR (200 perm)` = system.time(m_rwr(qq, ctx))[3],
      `Proximity-z (1000 perm)` = system.time(m_guney(qq, ctx))[3],
      `NetFunrankR-v1 statistic (1000 perm)` = if (n <= 500 || lib == "KEGG")
        system.time(m_v1(qq, ctx, B = 1000L))[3] else NA_real_)
    rt[[length(rt) + 1]] <- data.table(method = sub("\\.elapsed$", "", names(tm)),
                                       library = lib, size = n,
                                       n_perm = NA, seconds = unname(tm))
    print(rt[[length(rt)]])
  }
}
rt <- rbindlist(c(list(v1_time), rt), fill = TRUE)
fwrite(rt, file.path(RES, "runtime.csv"))

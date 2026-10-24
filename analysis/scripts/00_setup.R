# Shared setup for all benchmark scripts.
# Run from the package root:  Rscript analysis/scripts/<script>.R
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
})

AN <- "analysis"
# Private library for the version of the package under analysis, kept in the
# user cache directory (outside synced folders, whose file locks can break
# R CMD INSTALL). Override with the NF_LIB environment variable.
LIB <- Sys.getenv("NF_LIB",
                  file.path(tools::R_user_dir("NetFunrankR", "cache"), "analysis-lib"))
dir.create(LIB, recursive = TRUE, showWarnings = FALSE)

# Number of worker processes: NF_CORES, or all cores but one.
nf_cores <- function() {
  as.integer(Sys.getenv("NF_CORES", max(1L, parallel::detectCores() - 1L)))
}
.libPaths(c(normalizePath(LIB), .libPaths()))
if (!file.exists(file.path(LIB, "NetFunrankR", "DESCRIPTION")) ||
    isTRUE(as.logical(Sys.getenv("NF_REINSTALL", "FALSE")))) {
  system2(file.path(R.home("bin"), "R"),
          c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
            paste0("--library=", shQuote(normalizePath(LIB))), "."),
          stdout = FALSE)
}
suppressPackageStartupMessages(library(NetFunrankR, lib.loc = LIB))
stopifnot(exists("netEnrich", envir = asNamespace("NetFunrankR")))
data(COPD_proteinCodingGenes, package = "NetFunrankR")

# Thin wrappers returning data.frames, used throughout the benchmark scripts
nfTest <- function(query, gs, net, D = NULL, min_size = 5L, max_size = Inf,
                   background = NULL, beyond = "linker", exact_max = 500L) {
  as.data.frame(netEnrich(query, gs, net, distances = D, minSize = min_size,
                          maxSize = max_size, background = background,
                          beyond = beyond, exactMax = exact_max))
}
nfPerm <- function(query, gs, net, D = NULL, n_perm = 1000L,
                   statistic = "closest", min_size = 5L, max_size = Inf) {
  as.data.frame(netEnrichPerm(query, gs, net, distances = D, nPerm = n_perm,
                              statistic = statistic, minSize = min_size,
                              maxSize = max_size))
}

RES <- file.path(AN, "results")
FIG <- file.path(AN, "figures")
dir.create(RES, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)

# ---------------------------------------------------------------- networks
load_networks <- function() {
  nets <- readRDS(file.path(AN, "data", "networks.rds"))
  lapply(nets, buildNetwork)
}

# ---------------------------------------------------------------- gene sets
# KEGG (Enrichr 2021, symbols) for gene-list analyses
kegg_symbols <- function() loadPathways("KEGG")

# ---------------------------------------------------------------- competitors
# All competitors take (query, gs, net, ctx) where ctx caches per-library
# objects (incidence matrix, distance matrix) and return a named p-value
# vector over names(gs) (NA for untested sets).

make_ctx <- function(gs, net, min_size = 5L, max_size = 500L) {
  D <- geneSetDistances(net, gs)
  n_in <- colSums(D == 0L)
  keep <- names(n_in)[n_in >= min_size & n_in <= max_size]
  M <- Matrix::sparseMatrix(
    i = unlist(lapply(gs[keep], function(s) which(nodeNames(net) %in% s))),
    j = rep(seq_along(keep), vapply(gs[keep], function(s) sum(nodeNames(net) %in% s), 1L)),
    x = 1, dims = c(length(nodeNames(net)), length(keep)),
    dimnames = list(nodeNames(net), keep))
  W <- adjacency(net) %*% Matrix::Diagonal(x = 1 / nodeDegree(net))   # column-stochastic
  list(gs = gs, net = net, D = D, keep = keep, M = M, W = W,
       min_size = min_size, max_size = max_size,
       bins = split(seq_along(degreeBin(net)), degreeBin(net)))
}

# degree-bin preserving random replacements for query indices qi
draw_binned <- function(qi, ctx) {
  qb <- degreeBin(ctx$net)[qi]
  r <- integer(length(qi))
  for (b in unique(qb)) {
    w <- which(qb == b)
    pool <- ctx$bins[[as.character(b)]]
    r[w] <- pool[sample.int(length(pool), length(w))]
  }
  r
}

m_ora <- function(query, ctx) {
  # universe = all genes annotated in the library (standard ORA background)
  r <- oraTest(query, ctx$gs, minSize = 1L)
  stats::setNames(r$pvalue, r$set)[ctx$keep]
}

m_nfr <- function(query, ctx, background = NULL) {
  r <- nfTest(query, ctx$gs, ctx$net, D = ctx$D, background = background,
                         min_size = ctx$min_size, max_size = ctx$max_size)
  rd <- nfTest(query, ctx$gs, ctx$net, D = ctx$D, background = background,
                          min_size = ctx$min_size, max_size = ctx$max_size,
                          beyond = "distance")
  sv <- function(x, col) stats::setNames(x[[col]], x$set)
  list(`NetFunrankR-overlap` = sv(r, "p_overlap"), `NetFunrankR-linker` = sv(r, "p_linker"),
       `NetFunrankR-distance` = sv(r, "p_beyond"), `NetFunrankR-total` = sv(r, "p_proximity"),
       `NetFunrankR` = sv(r, "p_combined"), `NetFunrankR-distance-combined` = sv(rd, "p_combined"))
}

# NEA: edges between query and set vs configuration-model expectation
m_nea <- function(query, ctx) {
  net <- ctx$net
  x <- as.numeric(nodeNames(net) %in% query)
  nX <- as.numeric(adjacency(net) %*% x)                 # neighbours in X per node
  N <- as.numeric(Matrix::crossprod(ctx$M, nX))   # edges X-Y
  degX <- sum(nodeDegree(net) * x)
  degY <- as.numeric(Matrix::crossprod(ctx$M, nodeDegree(net)))
  E <- degX * degY / sum(nodeDegree(net))
  z <- (N - E) / sqrt(pmax(E, 1e-12))
  stats::setNames(stats::pnorm(z, lower.tail = FALSE), ctx$keep)
}

# RWR propagation; statistic = mean visiting probability over set members,
# null = degree-binned random seeds (B permutations, processed jointly)
m_rwr <- function(query, ctx, restart = 0.5, B = 200L, iters = 30L) {
  net <- ctx$net
  qi <- which(nodeNames(net) %in% query)
  seeds <- c(list(qi), lapply(seq_len(B), function(i) draw_binned(qi, ctx)))
  P0 <- Matrix::sparseMatrix(i = unlist(seeds),
                             j = rep(seq_along(seeds), lengths(seeds)),
                             x = 1 / length(qi),
                             dims = c(length(nodeNames(net)), length(seeds)))
  P0 <- as.matrix(P0)
  R <- P0
  for (it in seq_len(iters)) R <- (1 - restart) * as.matrix(ctx$W %*% R) + restart * P0
  S <- as.matrix(Matrix::crossprod(ctx$M, R)) / Matrix::colSums(ctx$M)
  obs <- S[, 1]
  p <- (1 + rowSums(S[, -1, drop = FALSE] >= obs)) / (B + 1)
  stats::setNames(p, ctx$keep)
}

# Closest-distance proximity d_c (Guney et al. 2016) with z-score and normal p
m_prox_closest <- function(query, ctx, B = 1000L) {
  r <- nfPerm(query, ctx$gs, ctx$net, D = ctx$D, n_perm = B,
                         statistic = "closest",
                         min_size = ctx$min_size, max_size = ctx$max_size)
  stats::setNames(r$p_norm, r$set)
}

# Symmetric module-to-module proximity d_XY (Menche et al. 2015) with a degree-binned null
m_prox_symmetric <- function(query, ctx, B = 200L) {
  r <- nfPerm(query, ctx$gs, ctx$net, D = ctx$D, n_perm = B,
                         statistic = "symmetric",
                         min_size = ctx$min_size, max_size = ctx$max_size)
  stats::setNames(r$p_norm, r$set)
}

# run a set of gene-list methods, returning long data.table(method, set, p)
run_genelist_methods <- function(query, ctx,
                                 methods = c("ora", "nea", "rwr", "closest", "symmetric", "nfr"),
                                 timing = FALSE) {
  out <- list()
  tm <- list()
  for (m in methods) {
    t0 <- proc.time()[["elapsed"]]
    p <- switch(m,
                ora = list(ORA = m_ora(query, ctx)),
                nea = list(NEA = m_nea(query, ctx)),
                rwr = list(RWR = m_rwr(query, ctx)),
                closest = list(`Proximity-closest` = m_prox_closest(query, ctx)),
                symmetric = list(`Proximity-symmetric` = m_prox_symmetric(query, ctx)),
                nfr = m_nfr(query, ctx))
    tm[[m]] <- proc.time()[["elapsed"]] - t0
    for (n in names(p))
      out[[n]] <- data.table(method = n, set = names(p[[n]]), p = unname(p[[n]]))
  }
  res <- rbindlist(out)
  if (timing) attr(res, "timing") <- unlist(tm)
  res
}

METHOD_LEVELS <- c("ORA", "NEA", "RWR", "Proximity-closest", "Proximity-symmetric",
                   "NetFunrankR-overlap", "NetFunrankR-linker", "NetFunrankR-distance", "NetFunrankR-total",
                   "NetFunrankR-distance-combined", "NetFunrankR")

# Re-run only the methods in `run` and replace their rows in `old` results
# (lets NetFunrankR be updated without recomputing slow competitors).
nf_methods_env <- function(default) {
  m <- Sys.getenv("NF_METHODS", "")
  if (nzchar(m)) strsplit(m, ",")[[1]] else default
}
merge_results <- function(old, new) {
  drop <- old$method %in% unique(new$method) |
    (grepl("^NetFunrankR", old$method) & any(grepl("^NetFunrankR", new$method)))
  rbind(old[!drop], new, fill = TRUE)
}
# Figures colour methods by family (3 categorical slots of the reference
# palette, validated for all pairs); individual methods are named on axes.
FAMILY_COLS <- c(`Annotation / expression` = "#2a78d6",
                 `Network (published)` = "#1baf7a",
                 `NetFunrankR` = "#eb6834")
method_family <- function(m) {
  ifelse(grepl("^NetFunrankR", m), "NetFunrankR",
         ifelse(m %in% c("NEA", "RWR", "Proximity-closest", "Proximity-symmetric"),
                "Network (published)", "Annotation / expression"))
}
# line charts with several individual methods: adjacent-validated order
LINE_COLS <- c(ORA = "#2a78d6", `NetFunrankR` = "#eb6834", `Proximity-closest` = "#1baf7a",
               RWR = "#eda100", `NetFunrankR-overlap` = "#e87ba4", `NetFunrankR-linker` = "#008300",
               NEA = "#6250d6", `NetFunrankR-distance` = "#e34948")

theme_nf <- function(base = 10) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

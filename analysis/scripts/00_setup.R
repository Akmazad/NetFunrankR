# Shared setup for all benchmark scripts.
# Run from the package root:  Rscript analysis/scripts/<script>.R
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
})

AN <- "analysis"
# private library outside OneDrive (sync locks break R CMD INSTALL)
LIB <- Sys.getenv("NF_LIB", file.path(Sys.getenv("LOCALAPPDATA"), "NetFunrankR_analysis_lib"))
dir.create(LIB, showWarnings = FALSE)
.libPaths(c(normalizePath(LIB), .libPaths()))
if (!file.exists(file.path(LIB, "NetFunrankR", "DESCRIPTION")) ||
    isTRUE(as.logical(Sys.getenv("NF_REINSTALL", "FALSE")))) {
  system2(file.path(R.home("bin"), "R"),
          c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
            paste0("--library=", shQuote(normalizePath(LIB))), "."),
          stdout = FALSE)
}
suppressPackageStartupMessages(library(NetFunrankR, lib.loc = LIB))
stopifnot(utils::packageVersion("NetFunrankR", lib.loc = LIB) >= "2.0.0")

RES <- file.path(AN, "results")
FIG <- file.path(AN, "figures")
dir.create(RES, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)

# ---------------------------------------------------------------- networks
load_networks <- function() {
  nets <- readRDS(file.path(AN, "data", "networks.rds"))
  lapply(nets, nf_network)
}

# ---------------------------------------------------------------- gene sets
# KEGG (Enrichr 2021, symbols) for gene-list analyses
kegg_symbols <- function() load_pathways("KEGG")

# ---------------------------------------------------------------- competitors
# All competitors take (query, gs, net, ctx) where ctx caches per-library
# objects (incidence matrix, distance matrix) and return a named p-value
# vector over names(gs) (NA for untested sets).

make_ctx <- function(gs, net, min_size = 5L, max_size = 500L) {
  D <- set_distances(net, gs)
  n_in <- colSums(D == 0L)
  keep <- names(n_in)[n_in >= min_size & n_in <= max_size]
  M <- Matrix::sparseMatrix(
    i = unlist(lapply(gs[keep], function(s) which(net$genes %in% s))),
    j = rep(seq_along(keep), vapply(gs[keep], function(s) sum(net$genes %in% s), 1L)),
    x = 1, dims = c(length(net$genes), length(keep)),
    dimnames = list(net$genes, keep))
  W <- net$adj %*% Matrix::Diagonal(x = 1 / net$degree)   # column-stochastic
  list(gs = gs, net = net, D = D, keep = keep, M = M, W = W,
       min_size = min_size, max_size = max_size,
       bins = split(seq_along(net$bin), net$bin))
}

# degree-bin preserving random replacements for query indices qi
draw_binned <- function(qi, ctx) {
  qb <- ctx$net$bin[qi]
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
  r <- nf_ora(query, ctx$gs, min_size = 1L)
  stats::setNames(r$pvalue, r$set)[ctx$keep]
}

m_nfr <- function(query, ctx, background = NULL) {
  r <- nf_proximity_test(query, ctx$gs, ctx$net, D = ctx$D, background = background,
                         min_size = ctx$min_size, max_size = ctx$max_size)
  rd <- nf_proximity_test(query, ctx$gs, ctx$net, D = ctx$D, background = background,
                          min_size = ctx$min_size, max_size = ctx$max_size,
                          beyond = "distance")
  sv <- function(x, col) stats::setNames(x[[col]], x$set)
  list(`NFR-overlap` = sv(r, "p_overlap"), `NFR-linker` = sv(r, "p_linker"),
       `NFR-distance` = sv(r, "p_beyond"), `NFR-total` = sv(r, "p_proximity"),
       `NFR-combined` = sv(r, "p_combined"), `NFR-combined-dist` = sv(rd, "p_combined"))
}

# NEA: edges between query and set vs configuration-model expectation
m_nea <- function(query, ctx) {
  net <- ctx$net
  x <- as.numeric(net$genes %in% query)
  nX <- as.numeric(net$adj %*% x)                 # neighbours in X per node
  N <- as.numeric(Matrix::crossprod(ctx$M, nX))   # edges X-Y
  degX <- sum(net$degree * x)
  degY <- as.numeric(Matrix::crossprod(ctx$M, net$degree))
  E <- degX * degY / sum(net$degree)
  z <- (N - E) / sqrt(pmax(E, 1e-12))
  stats::setNames(stats::pnorm(z, lower.tail = FALSE), ctx$keep)
}

# RWR propagation; statistic = mean visiting probability over set members,
# null = degree-binned random seeds (B permutations, processed jointly)
m_rwr <- function(query, ctx, restart = 0.5, B = 200L, iters = 30L) {
  net <- ctx$net
  qi <- which(net$genes %in% query)
  seeds <- c(list(qi), lapply(seq_len(B), function(i) draw_binned(qi, ctx)))
  P0 <- Matrix::sparseMatrix(i = unlist(seeds),
                             j = rep(seq_along(seeds), lengths(seeds)),
                             x = 1 / length(qi),
                             dims = c(length(net$genes), length(seeds)))
  P0 <- as.matrix(P0)
  R <- P0
  for (it in seq_len(iters)) R <- (1 - restart) * as.matrix(ctx$W %*% R) + restart * P0
  S <- as.matrix(Matrix::crossprod(ctx$M, R)) / Matrix::colSums(ctx$M)
  obs <- S[, 1]
  p <- (1 + rowSums(S[, -1, drop = FALSE] >= obs)) / (B + 1)
  stats::setNames(p, ctx$keep)
}

# Guney et al. (2016) closest-distance proximity with z-score and normal p
m_guney <- function(query, ctx, B = 1000L) {
  r <- nf_proximity_perm(query, ctx$gs, ctx$net, D = ctx$D, n_perm = B,
                         statistic = "closest",
                         min_size = ctx$min_size, max_size = ctx$max_size)
  stats::setNames(r$p_norm, r$set)
}

# NetFunrankR v1 statistic (symmetric d_XY) with corrected binned null
m_v1 <- function(query, ctx, B = 200L) {
  r <- nf_proximity_perm(query, ctx$gs, ctx$net, D = ctx$D, n_perm = B,
                         statistic = "symmetric",
                         min_size = ctx$min_size, max_size = ctx$max_size)
  stats::setNames(r$p_norm, r$set)
}

# run a set of gene-list methods, returning long data.table(method, set, p)
run_genelist_methods <- function(query, ctx,
                                 methods = c("ora", "nea", "rwr", "guney", "v1", "nfr"),
                                 timing = FALSE) {
  out <- list()
  tm <- list()
  for (m in methods) {
    t0 <- proc.time()[["elapsed"]]
    p <- switch(m,
                ora = list(ORA = m_ora(query, ctx)),
                nea = list(NEA = m_nea(query, ctx)),
                rwr = list(RWR = m_rwr(query, ctx)),
                guney = list(`Proximity-z` = m_guney(query, ctx)),
                v1 = list(`NetFunrankR-v1` = m_v1(query, ctx)),
                nfr = m_nfr(query, ctx))
    tm[[m]] <- proc.time()[["elapsed"]] - t0
    for (n in names(p))
      out[[n]] <- data.table(method = n, set = names(p[[n]]), p = unname(p[[n]]))
  }
  res <- rbindlist(out)
  if (timing) attr(res, "timing") <- unlist(tm)
  res
}

METHOD_LEVELS <- c("ORA", "NEA", "RWR", "Proximity-z", "NetFunrankR-v1",
                   "NFR-overlap", "NFR-linker", "NFR-distance", "NFR-total",
                   "NFR-combined-dist", "NFR-combined")

# Re-run only the methods in `run` and replace their rows in `old` results
# (lets NetFunrankR be updated without recomputing slow competitors).
nf_methods_env <- function(default) {
  m <- Sys.getenv("NF_METHODS", "")
  if (nzchar(m)) strsplit(m, ",")[[1]] else default
}
merge_results <- function(old, new) {
  drop <- old$method %in% unique(new$method) |
    (grepl("^NFR", old$method) & any(grepl("^NFR", new$method)))
  rbind(old[!drop], new, fill = TRUE)
}
# Figures colour methods by family (3 categorical slots of the reference
# palette, validated for all pairs); individual methods are named on axes.
FAMILY_COLS <- c(`Annotation / expression` = "#2a78d6",
                 `Network (published)` = "#1baf7a",
                 `NetFunrankR v2` = "#eb6834")
method_family <- function(m) {
  ifelse(grepl("^NFR", m), "NetFunrankR v2",
         ifelse(m %in% c("NEA", "RWR", "Proximity-z", "NetFunrankR-v1"),
                "Network (published)", "Annotation / expression"))
}
# line charts with several individual methods: adjacent-validated order
LINE_COLS <- c(ORA = "#2a78d6", `NFR-combined` = "#eb6834", `Proximity-z` = "#1baf7a",
               RWR = "#eda100", `NFR-overlap` = "#e87ba4", `NFR-linker` = "#008300",
               NEA = "#6250d6", `NFR-distance` = "#e34948")

theme_nf <- function(base = 10) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

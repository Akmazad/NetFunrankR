# Decomposed network-proximity enrichment test (exact, permutation-free)
#
# Notation (see manuscript, Methods):
#   X  query genes in the network,   Y  gene set,   d_Y(v) = min_y d(v, y)
#   S  = sum_{x in X} d_Y(x)                  (closest-distance proximity, d_c)
#   O  = |X intersect Y|                       (overlap)
#   S' = sum_{x in X \ Y} d_Y(x)               (beyond-overlap proximity; S' = S)
# Null: every query gene is replaced by a node drawn uniformly from its own
# degree bin. S, O are then sums of independent bin-specific variables whose
# exact distributions are obtained by direct discrete convolution
# (src/exact_null.cpp). All terms are non-negative, so tails down to ~1e-300
# are accurate, unlike FFT-based convolution.

# Distribution summaries for a pmf on 0..(length-1)
.lower_p <- function(p, s) {          # P(S <= s), P(S < s) + P(S = s)/2
  cp <- cumsum(p)
  i <- s + 1L
  c(p = min(1, cp[i]), mid = min(1, cp[i] - p[i] / 2))
}
.upper_p <- function(p, o) {          # P(O >= o), P(O > o) + P(O = o)/2
  tail <- sum(p[(o + 1L):length(p)])
  c(p = min(1, tail), mid = min(1, tail - p[o + 1L] / 2))
}

#' Decomposed network-proximity enrichment test
#'
#' For every gene set \eqn{Y}, tests whether the query genes \eqn{X} are
#' closer to \eqn{Y} in the interactome than degree-matched random genes, and
#' decomposes this proximity into (i) direct membership (overlap) and
#' (ii) network proximity of the non-member query genes ("beyond overlap").
#' All null distributions are exact (discrete convolution over degree bins);
#' no permutations are required.
#'
#' Reported p-values (one-sided):
#' \describe{
#'   \item{p_overlap}{degree-aware overlap test: \eqn{P(O \ge o)} with
#'     \eqn{O} a Poisson-binomial sum over query genes of the probability that
#'     a node of the same degree bin belongs to \eqn{Y}.}
#'   \item{p_linker}{conditional linker test: \eqn{P(L \ge l)} for the number
#'     \eqn{L} of non-member query genes adjacent to \eqn{Y}, where each
#'     non-member is replaced by a random non-member of its degree bin
#'     (Poisson-binomial). Uniform given the membership pattern, hence
#'     independent of \code{p_overlap} under the null.}
#'   \item{p_beyond}{conditional distance test of non-member query genes:
#'     \eqn{P(S' \le s')} under the same conditional null.}
#'   \item{p_proximity}{total closest-distance proximity \eqn{P(S \le s)}.}
#'   \item{p_combined}{Fisher combination of the mid-p values of
#'     \code{p_overlap} and the beyond-overlap test selected by \code{beyond}
#'     (\code{p_linker} by default; chi-square, 4 df).}
#'   \item{p_hyper}{classical hypergeometric ORA, for reference.}
#' }
#'
#' @param query character vector of query genes (same ID type as the network).
#' @param genesets named list of gene vectors, or an annotation data.frame.
#' @param net an \code{nf_network} (see \code{\link{nf_network}}).
#' @param D optional precomputed \code{\link{set_distances}} matrix for
#'   \code{net} and \code{genesets} (re-use it when testing many queries).
#' @param min_size,max_size gene sets with fewer/more genes in the network are
#'   skipped.
#' @param universe background for the hypergeometric ORA; default all genes
#'   annotated in \code{genesets}.
#' @param background optional set of genes that could have been selected
#'   (e.g. the genes measured in an expression experiment). Null replacements
#'   are then drawn from background genes of the same degree bin only.
#' @param beyond beyond-overlap test entering \code{p_combined}: "linker"
#'   (count of adjacent non-member query genes; robust when only a few query
#'   genes relate to the set) or "distance" (sum of distances).
#' @param exact_max query lists with at most this many genes in the network
#'   use exact convolution; larger lists use a continuity-corrected
#'   saddlepoint approximation (Daniels 1987) of the same null distributions.
#' @param p_adjust multiple-testing method passed to \code{p.adjust}.
#'
#' @returns data.frame with one row per tested gene set, sorted by
#'   \code{p_combined}.
#' @export
nf_proximity_test <- function(query, genesets, net, D = NULL,
                              min_size = 5L, max_size = Inf,
                              universe = NULL, background = NULL,
                              beyond = c("linker", "distance"),
                              exact_max = 500L, p_adjust = "BH") {
  beyond <- match.arg(beyond)
  genesets <- as_geneset_list(genesets)
  query <- unique(as.character(query))
  if (is.null(D)) D <- set_distances(net, genesets)
  stopifnot(identical(rownames(D), net$genes))
  genesets <- genesets[colnames(D)]

  qi <- which(net$genes %in% query)
  if (length(qi) < 2L)
    stop("fewer than 2 query genes are in the network")
  qbin <- net$bin[qi]
  B <- max(net$bin)
  m_b <- tabulate(qbin, nbins = B)

  if (is.null(universe)) universe <- unique(unlist(genesets, use.names = FALSE))
  q_univ <- intersect(query, universe)

  # nodes that random replacements may be drawn from
  bg <- if (is.null(background)) seq_along(net$genes) else
    sort(union(which(net$genes %in% background), qi))
  bg_bin <- net$bin[bg]

  n_in <- colSums(D == 0L, na.rm = TRUE)
  keep <- which(n_in >= min_size & n_in <= max_size)

  res <- lapply(keep, function(k) {
    d <- D[, k]
    dmax <- max(d)
    # bin x distance counts over the background, then row-normalised kernels
    cnt <- matrix(tabulate(d[bg] * B + bg_bin, nbins = B * (dmax + 1L)),
                  nrow = B, ncol = dmax + 1L)
    f <- cnt / pmax(rowSums(cnt), 1)

    dq <- d[qi]
    member <- dq == 0L
    o <- sum(member)
    s <- sum(dq)
    mnon_b <- tabulate(qbin[!member], nbins = B)

    # kernels: overlap (Bernoulli), total distance (f), conditional distance
    # of non-members (g = f restricted to d >= 1)
    fO <- cbind(1 - f[, 1L], f[, 1L])
    g <- f
    g[, 1L] <- 0
    gs_row <- rowSums(g)
    g <- g / pmax(gs_row, 1e-300)
    g[gs_row == 0, 2L] <- 1

    # linker count: non-member query genes adjacent to Y; conditional on the
    # membership pattern each non-member is a linker with probability g[b, 2]
    lnk <- sum(dq == 1L)
    fL <- cbind(1 - g[, 2L], g[, 2L])

    if (length(qi) <= exact_max) {
      nulls <- exact_null_cpp(f, m_b, mnon_b)
      po <- .upper_p(nulls$O, o)
      ps <- .lower_p(nulls$S, s)
      pb <- if (o < length(qi)) .lower_p(nulls$Sb, s) else c(p = 1, mid = 1)
      pl <- if (o < length(qi))
        # exact_null_cpp()$O counts successes with probability in column 1
        .upper_p(exact_null_cpp(fL[, 2:1, drop = FALSE], mnon_b, integer(B))$O, lnk)
      else c(p = 1, mid = 1)
    } else {
      sp_up <- function(k) .sp_upper(fO, m_b, k)
      po <- c(p = sp_up(o), mid = (sp_up(o) + sp_up(o + 1L)) / 2)
      sp_lo <- function(f, m, x) .sp_lower(f, m, x)
      ps <- c(p = sp_lo(f, m_b, s), mid = (sp_lo(f, m_b, s) + sp_lo(f, m_b, s - 1L)) / 2)
      pb <- if (o < length(qi))
        c(p = sp_lo(g, mnon_b, s), mid = (sp_lo(g, mnon_b, s) + sp_lo(g, mnon_b, s - 1L)) / 2)
      else c(p = 1, mid = 1)
      pl <- if (o < length(qi)) {
        a <- .sp_upper(fL, mnon_b, lnk)
        c(p = a, mid = (a + .sp_upper(fL, mnon_b, lnk + 1L)) / 2)
      } else c(p = 1, mid = 1)
    }
    mom_S <- .sum_moments(f, m_b)
    if (o < length(qi)) {
      mom_b <- .sum_moments(g, mnon_b)
      e_b <- mom_b[["mean"]] / (length(qi) - o)
      z_b <- (mom_b[["mean"]] - s) / max(mom_b[["sd"]], 1e-12)
    } else {
      e_b <- NA_real_
      z_b <- NA_real_
    }
    pc <- if (beyond == "linker") pl else pb
    chi <- -2 * (log(max(po[["mid"]], 1e-300)) + log(max(pc[["mid"]], 1e-300)))

    # classical hypergeometric ORA
    Yset <- genesets[[k]]
    M <- length(intersect(Yset, universe))
    N <- length(universe)
    xh <- length(intersect(q_univ, Yset))
    p_h <- stats::phyper(xh - 1L, M, N - M, length(q_univ), lower.tail = FALSE)

    data.frame(
      set = colnames(D)[k],
      set_size = length(Yset),
      set_in_net = n_in[[k]],
      n_query_net = length(qi),
      overlap = o,
      overlap_expected = sum(m_b * f[, 1L]),
      d_c = s / length(qi),
      d_c_expected = mom_S[["mean"]] / length(qi),
      z_proximity = (mom_S[["mean"]] - s) / max(mom_S[["sd"]], 1e-12),
      d_beyond = if (o < length(qi)) s / (length(qi) - o) else NA_real_,
      d_beyond_expected = e_b,
      z_beyond = z_b,
      n_linker = lnk,
      linker_expected = sum(mnon_b * g[, 2L]),
      p_hyper = p_h,
      p_overlap = po[["p"]],
      p_linker = pl[["p"]],
      p_beyond = pb[["p"]],
      p_proximity = ps[["p"]],
      p_combined = stats::pchisq(chi, df = 4, lower.tail = FALSE),
      overlap_genes = paste(net$genes[qi][member], collapse = "/"),
      linker_genes = paste(net$genes[qi][dq == 1L], collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, res)
  for (p in c("p_hyper", "p_overlap", "p_linker", "p_beyond", "p_proximity", "p_combined"))
    res[[sub("^p_", "fdr_", p)]] <- stats::p.adjust(res[[p]], method = p_adjust)
  res <- res[order(res$p_combined, res$p_proximity), ]
  rownames(res) <- NULL
  attr(res, "query_not_in_network") <- setdiff(query, net$genes)
  class(res) <- c("nf_result", class(res))
  res
}

#' Classify gene sets by the source of their enrichment evidence
#'
#' Sets significant by the combined test (FDR < \code{alpha}) are labelled by
#' which component reaches nominal significance (p < \code{alpha}).
#'
#' @param res an \code{nf_result} from \code{\link{nf_proximity_test}}.
#' @param alpha significance threshold.
#' @returns factor: "overlap + network", "overlap only", "network only",
#'   "joint" (only the combination is significant) or "none".
#' @export
nf_evidence_class <- function(res, alpha = 0.05) {
  o <- res$p_overlap < alpha
  b <- res$p_linker < alpha
  out <- ifelse(o & b, "overlap + network",
                ifelse(o, "overlap only", ifelse(b, "network only", "joint")))
  out[res$fdr_combined >= alpha] <- "none"
  factor(out, levels = c("overlap + network", "overlap only", "network only",
                         "joint", "none"))
}

#' Monte-Carlo permutation version of the proximity test
#'
#' Degree-bin-preserving permutation null, used (i) to validate the exact
#' convolution null and (ii) to compute the symmetric module-to-module
#' proximity of NetFunrankR v1 (Eq. 2 of Menche et al. 2015 / Guney et al.
#' 2016) for comparison.
#'
#' @inheritParams nf_proximity_test
#' @param n_perm number of permutations.
#' @param statistic "closest" (\eqn{d_c}) or "symmetric" (\eqn{d_{XY}}).
#' @param batch permutations processed together in one multi-source BFS
#'   (symmetric statistic only).
#' @returns data.frame with observed statistic, null mean/sd, z-score,
#'   empirical p-value \eqn{(1 + \#\{null \le obs\})/(1 + B)} and the
#'   normal-approximation p-value used in v1.
#' @export
nf_proximity_perm <- function(query, genesets, net, D = NULL, n_perm = 1000L,
                              statistic = c("closest", "symmetric"),
                              min_size = 5L, max_size = Inf, batch = 50L) {
  statistic <- match.arg(statistic)
  genesets <- as_geneset_list(genesets)
  if (is.null(D)) D <- set_distances(net, genesets)
  n_in <- colSums(D == 0L, na.rm = TRUE)
  keep <- which(n_in >= min_size & n_in <= max_size)
  D <- D[, keep, drop = FALSE]
  qi <- which(net$genes %in% query)
  bins <- split(seq_along(net$bin), net$bin)
  qbin <- net$bin[qi]

  draw <- function() {
    r <- integer(length(qi))
    for (b in unique(qbin)) {
      w <- which(qbin == b)
      pool <- bins[[as.character(b)]]
      r[w] <- pool[sample.int(length(pool), length(w))]
    }
    r
  }

  member <- D == 0L
  ysz <- colSums(member)
  sym_stat <- function(xi, dX) {
    # dX: distances from X to every node (vector)
    (colSums(D[xi, , drop = FALSE]) + colSums(member * dX)) / (length(xi) + ysz)
  }

  if (statistic == "closest") {
    obs <- colMeans(D[qi, , drop = FALSE])
    null <- vapply(seq_len(n_perm), function(i) colMeans(D[draw(), , drop = FALSE]),
                   numeric(ncol(D)))
  } else {
    dq <- set_distances(net, list(q = net$genes[qi]))[, 1]
    obs <- sym_stat(qi, dq)
    null <- matrix(NA_real_, ncol(D), n_perm)
    for (st in seq(1L, n_perm, by = batch)) {
      ids <- st:min(n_perm, st + batch - 1L)
      R <- lapply(ids, function(i) draw())
      dR <- set_distances(net, stats::setNames(lapply(R, function(r) net$genes[r]),
                                               ids))
      for (j in seq_along(ids)) null[, ids[j]] <- sym_stat(R[[j]], dR[, j])
    }
  }
  if (!is.matrix(null)) null <- matrix(null, nrow = ncol(D))
  mu <- rowMeans(null)
  sdv <- apply(null, 1, stats::sd)
  z <- (obs - mu) / sdv
  data.frame(set = colnames(D), statistic = statistic, observed = obs,
             null_mean = mu, null_sd = sdv, z = z,
             p_perm = (1 + rowSums(null <= obs)) / (1 + n_perm),
             p_norm = stats::pnorm(z),
             stringsAsFactors = FALSE)
}

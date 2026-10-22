# Saddlepoint tail probabilities for sums of independent integer-valued
# variables (used for large query lists, where exact convolution is costly).
#
# S = sum over bins b of m_b i.i.d. draws from the kernel f[b, ] on 0..L-1.
# Upper tail P(S >= k) by the second continuity-corrected Lugannani-Rice
# formula for lattice variables (Daniels 1987, Int Stat Rev 55:37-48):
#   K'(t) = k - 1/2,  w = sgn(t) sqrt(2 (t (k - 1/2) - K(t))),
#   u = 2 sinh(t/2) sqrt(K''(t)),  P ~ 1 - Phi(w) - phi(w) (1/w - 1/u).

.sp_cgf <- function(t, f, m) {
  j <- seq_len(ncol(f)) - 1
  lw <- sweep(log(f), 2, t * j, "+")              # log f_bj + t j
  mx <- apply(lw, 1, max)
  w <- exp(lw - mx)
  z <- rowSums(w)
  p <- w / z                                       # tilted kernels
  e1 <- as.numeric(p %*% j)
  e2 <- as.numeric(p %*% j^2)
  c(K = sum(m * (mx + log(z))), K1 = sum(m * e1), K2 = sum(m * (e2 - e1^2)))
}

# P(S >= k)
.sp_upper <- function(f, m, k) {
  keep <- m > 0
  f <- f[keep, , drop = FALSE]
  m <- m[keep]
  supp <- f > 0
  lo <- sum(m * apply(supp, 1, function(r) min(which(r)) - 1))
  hi <- sum(m * apply(supp, 1, function(r) max(which(r)) - 1))
  if (k <= lo) return(1)
  if (k > hi) return(0)
  if (k == hi) {
    top <- apply(supp, 1, function(r) max(which(r)))
    return(exp(sum(m * log(f[cbind(seq_len(nrow(f)), top)]))))
  }
  kt <- k - 0.5
  # Newton iterations safeguarded by bisection on the monotone K'(t)
  a <- -60
  b <- 60
  t <- 0
  for (it in seq_len(200)) {
    cg <- .sp_cgf(t, f, m)
    g <- cg[["K1"]] - kt
    if (abs(g) < 1e-10 * max(1, abs(kt))) break
    if (g > 0) b <- t else a <- t
    tn <- t - g / cg[["K2"]]
    t <- if (is.finite(tn) && tn > a && tn < b) tn else (a + b) / 2
  }
  cg <- .sp_cgf(t, f, m)
  if (abs(t) < 1e-5) {
    mu <- .sp_cgf(0, f, m)
    return(stats::pnorm((kt - mu[["K1"]]) / sqrt(mu[["K2"]]), lower.tail = FALSE))
  }
  w <- sign(t) * sqrt(max(2 * (t * kt - cg[["K"]]), 0))
  u <- 2 * sinh(t / 2) * sqrt(cg[["K2"]])
  p <- stats::pnorm(w, lower.tail = FALSE) - stats::dnorm(w) * (1 / w - 1 / u)
  min(max(p, 0), 1)
}

# P(S <= s) via the reflected kernel
.sp_lower <- function(f, m, s) {
  L <- ncol(f)
  .sp_upper(f[, L:1, drop = FALSE], m, sum(m) * (L - 1) - s)
}

# moments of S without forming its distribution
.sum_moments <- function(f, m) {
  j <- seq_len(ncol(f)) - 1
  e1 <- as.numeric(f %*% j)
  e2 <- as.numeric(f %*% j^2)
  c(mean = sum(m * e1), sd = sqrt(sum(m * (e2 - e1^2))))
}

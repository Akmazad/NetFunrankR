toy_network <- function(n = 300, m = 3, seed = 1) {
  set.seed(seed)
  g <- igraph::sample_pa(n, m = m, directed = FALSE)
  igraph::V(g)$name <- paste0("G", seq_len(n))
  nf_network(igraph::as_data_frame(g)[, 1:2], min_bin_size = 20L)
}

test_that("degree bins are monotone and respect the minimum size", {
  deg <- c(rep(1L, 50), rep(2L, 30), rep(3L, 10), rep(10L, 5), 50L)
  b <- degree_bins(deg, min_bin_size = 20L)
  expect_true(all(diff(b[order(deg)]) >= 0))
  expect_true(all(table(b) >= 20L))
  expect_equal(length(unique(b[deg == 2L])), 1L)
})

test_that("set_distances matches igraph shortest paths", {
  net <- toy_network()
  gs <- list(A = paste0("G", 1:5), B = paste0("G", 100:120))
  D <- set_distances(net, gs)
  ref <- apply(igraph::distances(net$graph, v = gs$B, to = net$genes), 2, min)
  expect_equal(unname(D[, "B"]), unname(as.integer(ref)))
})

test_that("ORA equals the one-sided Fisher exact test", {
  universe <- paste0("G", 1:200)
  gs <- list(S = paste0("G", 1:30))
  q <- paste0("G", c(1:6, 50:60))
  r <- nf_ora(q, gs, universe = universe, min_size = 1L)
  tab <- matrix(c(6, 24, 11, 159), 2)
  expect_equal(r$pvalue, stats::fisher.test(tab, alternative = "greater")$p.value)
})

test_that("exact null equals brute-force enumeration", {
  f <- rbind(c(0.2, 0.5, 0.3), c(0.1, 0.1, 0.8))
  m <- c(2L, 1L)
  nl <- exact_null_cpp(f, m, c(0L, 0L))
  grid <- expand.grid(a = 0:2, b = 0:2, c = 0:2)
  pr <- f[1, grid$a + 1] * f[1, grid$b + 1] * f[2, grid$c + 1]
  bf <- tapply(pr, grid$a + grid$b + grid$c, sum)
  expect_equal(nl$S, as.numeric(bf))
  expect_equal(sum(nl$O), 1)
  expect_equal(nl$O, c(0.8^2 * 0.9, 2 * 0.2 * 0.8 * 0.9 + 0.8^2 * 0.1,
                       0.2^2 * 0.9 + 2 * 0.2 * 0.8 * 0.1, 0.2^2 * 0.1))
})

test_that("exact test agrees with Monte-Carlo permutation", {
  net <- toy_network()
  set.seed(2)
  gs <- list(A = sample(net$genes, 25), B = sample(net$genes, 40))
  q <- c(gs$A[1:3], sample(net$genes, 20))
  ex <- nf_proximity_test(q, gs, net, min_size = 1L)
  mc <- nf_proximity_perm(q, gs, net, n_perm = 4000L, min_size = 1L)
  m <- merge(ex, mc, by = "set")
  expect_true(all(abs(m$p_proximity - m$p_perm) < 0.03))
})

test_that("saddlepoint approximates the exact null", {
  net <- toy_network(n = 600)
  set.seed(3)
  gs <- list(A = sample(net$genes, 40), B = sample(net$genes, 60))
  q <- sample(net$genes, 120)
  ex <- nf_proximity_test(q, gs, net, min_size = 1L, exact_max = 1e6)
  sp <- nf_proximity_test(q, gs, net, min_size = 1L, exact_max = 0)
  m <- merge(ex, sp, by = "set")
  for (p in c("p_overlap", "p_beyond", "p_proximity"))
    expect_true(all(abs(log10(m[[paste0(p, ".x")]]) - log10(m[[paste0(p, ".y")]])) < 0.05))
})

test_that("beyond-overlap p-values are uniform under the degree-matched null", {
  net <- toy_network(n = 500)
  set.seed(4)
  gs <- list(A = sample(net$genes, 50))
  p <- replicate(300, nf_proximity_test(sample(net$genes, 20), gs, net,
                                        min_size = 1L)$p_beyond)
  expect_lt(mean(p < 0.05), 0.08)
})

test_that("linker test is calibrated, detects planted linkers and matches saddlepoint", {
  net <- toy_network(n = 600)
  set.seed(5)
  gs <- list(A = sample(net$genes, 40))
  p <- replicate(300, nf_proximity_test(sample(net$genes, 25), gs, net,
                                        min_size = 1L)$p_linker)
  expect_lt(mean(p < 0.05), 0.08)
  net2 <- toy_network(n = 1000, m = 2, seed = 6)
  gs2 <- list(A = sample(net2$genes, 10))
  nb <- unique(unlist(igraph::adjacent_vertices(net2$graph, gs2$A)))
  linkers <- setdiff(igraph::V(net2$graph)$name[nb], gs2$A)
  q <- c(sample(linkers, 8), sample(setdiff(net2$genes, c(gs2$A, linkers)), 17))
  r <- nf_proximity_test(q, gs2, net2, min_size = 1L)
  expect_lt(r$p_linker, 1e-3)
  sp <- nf_proximity_test(q, gs2, net2, min_size = 1L, exact_max = 0)
  expect_lt(abs(log10(r$p_linker) - log10(sp$p_linker)), 0.1)
})

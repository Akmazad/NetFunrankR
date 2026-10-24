test_that("degree bins are monotone and respect the minimum size", {
    deg <- c(rep(1L, 50), rep(2L, 30), rep(3L, 10), rep(10L, 5), 50L)
    b <- degreeBins(deg, minBinSize = 20L)
    expect_true(all(diff(b[order(deg)]) >= 0))
    expect_true(all(table(b) >= 20L))
    expect_length(unique(b[deg == 2L]), 1L)
})

test_that("geneSetDistances matches igraph shortest paths", {
    net <- toyNetwork()
    gs <- list(A = paste0("G", 1:5), B = paste0("G", 100:120))
    D <- geneSetDistances(net, gs)
    ref <- apply(igraph::distances(nfGraph(net), v = gs$B,
        to = nodeNames(net)), 2, min)
    expect_equal(unname(D[, "B"]), unname(as.integer(ref)))
})

test_that("oraTest equals the one-sided Fisher exact test", {
    universe <- paste0("G", 1:200)
    gs <- list(S = paste0("G", 1:30))
    q <- paste0("G", c(1:6, 50:60))
    r <- oraTest(q, gs, universe = universe, minSize = 1L)
    tab <- matrix(c(6, 24, 11, 159), 2)
    expect_equal(r$pvalue,
        stats::fisher.test(tab, alternative = "greater")$p.value)
})

test_that("exact null equals brute-force enumeration", {
    f <- rbind(c(0.2, 0.5, 0.3), c(0.1, 0.1, 0.8))
    nl <- exact_null_cpp(f, c(2L, 1L), c(0L, 0L))
    grid <- expand.grid(a = 0:2, b = 0:2, c = 0:2)
    pr <- f[1, grid$a + 1] * f[1, grid$b + 1] * f[2, grid$c + 1]
    expect_equal(nl$S, as.numeric(tapply(pr, grid$a + grid$b + grid$c, sum)))
    expect_equal(nl$O, c(0.8^2 * 0.9, 2 * 0.2 * 0.8 * 0.9 + 0.8^2 * 0.1,
        0.2^2 * 0.9 + 2 * 0.2 * 0.8 * 0.1, 0.2^2 * 0.1))
})

test_that("exact test agrees with Monte-Carlo permutation", {
    net <- toyNetwork()
    set.seed(2)
    gs <- list(A = sample(nodeNames(net), 25), B = sample(nodeNames(net), 40))
    q <- c(gs$A[1:3], sample(nodeNames(net), 20))
    ex <- as.data.frame(netEnrich(q, gs, net, minSize = 1L))
    mc <- as.data.frame(netEnrichPerm(q, gs, net, nPerm = 4000L,
        minSize = 1L))
    m <- merge(ex, mc, by = "set")
    expect_true(all(abs(m$p_proximity - m$p_perm) < 0.03))
})

test_that("saddlepoint approximates the exact null", {
    net <- toyNetwork(n = 600)
    set.seed(3)
    gs <- list(A = sample(nodeNames(net), 40), B = sample(nodeNames(net), 60))
    q <- sample(nodeNames(net), 120)
    ex <- as.data.frame(netEnrich(q, gs, net, minSize = 1L, exactMax = 1e6))
    sp <- as.data.frame(netEnrich(q, gs, net, minSize = 1L, exactMax = 0))
    m <- merge(ex, sp, by = "set")
    for (p in c("p_overlap", "p_linker", "p_beyond", "p_proximity")) {
        d <- abs(log10(m[[paste0(p, ".x")]]) - log10(m[[paste0(p, ".y")]]))
        expect_true(all(d < 0.05), label = p)
    }
})

test_that("beyond-overlap tests are calibrated under the null", {
    net <- toyNetwork(n = 600)
    set.seed(4)
    gs <- list(A = sample(nodeNames(net), 40))
    p <- replicate(300, unlist(as.data.frame(netEnrich(
        sample(nodeNames(net), 25), gs, net,
        minSize = 1L))[c("p_linker", "p_beyond", "p_combined")]))
    expect_true(all(rowMeans(p < 0.05) < 0.08))
})

test_that("planted linker genes are detected", {
    net <- toyNetwork(n = 1000, m = 2, seed = 6)
    gs <- list(A = sample(nodeNames(net), 10))
    nb <- unique(unlist(igraph::adjacent_vertices(nfGraph(net), gs$A)))
    linkers <- setdiff(nodeNames(net)[nb], gs$A)
    q <- c(sample(linkers, 8),
        sample(setdiff(nodeNames(net), c(gs$A, linkers)), 17))
    r <- enrichResults(netEnrich(q, gs, net, minSize = 1L))
    expect_lt(r$p_linker, 1e-3)
    expect_equal(r$n_linker, 8L)
    expect_identical(as.character(evidenceClass(netEnrich(q, gs, net,
        minSize = 1L))), "network only")
})

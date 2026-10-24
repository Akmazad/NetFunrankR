test_that("NFNetwork objects are valid and have working accessors", {
    net <- toyNetwork()
    expect_s4_class(net, "NFNetwork")
    expect_true(methods::validObject(net))
    expect_length(nodeNames(net), length(nodeDegree(net)))
    expect_equal(as.integer(Matrix::rowSums(adjacency(net))),
        unname(nodeDegree(net)))
    expect_output(show(net), "NFNetwork")
    expect_error(methods::new("NFNetwork", graph = 1), "igraph")
})

test_that("buildNetwork keeps the largest connected component", {
    edges <- data.frame(from = c("A", "B", "X"), to = c("B", "C", "Y"))
    net <- buildNetwork(edges, minBinSize = 1L)
    expect_setequal(nodeNames(net), c("A", "B", "C"))
    expect_error(buildNetwork(data.frame(a = "A", b = "A")), "no valid")
})

test_that("NFResult accessors, coercion and export work", {
    net <- toyNetwork()
    set.seed(7)
    gs <- list(A = sample(nodeNames(net), 30), B = sample(nodeNames(net), 30))
    q <- c(sample(nodeNames(net), 20), "NOT_A_GENE")
    res <- netEnrich(q, gs, net, minSize = 1L)
    expect_s4_class(res, "NFResult")
    expect_s4_class(enrichResults(res), "DataFrame")
    expect_identical(unmatchedGenes(res), "NOT_A_GENE")
    expect_length(queryGenes(res), 20L)
    expect_setequal(names(testedGeneSets(res)), c("A", "B"))
    df <- as.data.frame(res)
    expect_true(all(c("evidence", "p_combined", "fdr_combined") %in%
        colnames(df)))
    expect_true(all(levels(evidenceClass(res)) %in% names(evidenceColors())))
    f <- exportResults(res, tempfile(fileext = ".tsv"))
    expect_equal(nrow(utils::read.delim(f)), 2L)
    expect_output(show(res), "NFResult")
})

test_that("gene sets are accepted in all supported formats", {
    gs <- list(S1 = c("G1", "G2"), S2 = c("G3", NA, "G3"))
    expect_identical(asGeneSetList(gs)$S2, "G3")
    df <- data.frame(GeneSymbol = c("G1", "G2", "G3"),
        PathwayName = c("S1", "S1", "S2"))
    expect_identical(asGeneSetList(df), list(S1 = c("G1", "G2"), S2 = "G3"))
    gmt <- tempfile(fileext = ".gmt")
    writeLines(c("S1\tdesc\tG1\tG2", "S2\tdesc\tG3"), gmt)
    expect_identical(asGeneSetList(gmt), list(S1 = c("G1", "G2"), S2 = "G3"))
    skip_if_not_installed("GSEABase")
    gsc <- GSEABase::GeneSetCollection(list(
        GSEABase::GeneSet(c("G1", "G2"), setName = "S1")))
    expect_identical(asGeneSetList(gsc), list(S1 = c("G1", "G2")))
    expect_error(asGeneSetList(1:3), "named list")
})

test_that("netEnrichBatch matches single runs", {
    net <- toyNetwork()
    set.seed(8)
    gs <- list(A = sample(nodeNames(net), 30))
    q <- list(x = sample(nodeNames(net), 15), y = sample(nodeNames(net), 25))
    b <- netEnrichBatch(q, gs, net, minSize = 1L)
    expect_named(b, c("x", "y"))
    expect_equal(enrichResults(b$y)$p_combined,
        enrichResults(netEnrich(q$y, gs, net, minSize = 1L))$p_combined)
})

test_that("plots return ggplot objects", {
    net <- toyNetwork()
    set.seed(9)
    gs <- list(A = sample(nodeNames(net), 30), B = sample(nodeNames(net), 30))
    q <- c(gs$A[1:5], sample(nodeNames(net), 15))
    res <- netEnrich(q, gs, net, minSize = 1L)
    expect_s3_class(plotEvidenceMap(res), "ggplot")
    expect_s3_class(plotTopSets(res), "ggplot")
    expect_s3_class(plotLinkerNetwork(res, "A", net), "ggplot")
    nd <- nullDistribution(q, gs, net, "A")
    expect_equal(sum(nd$S), 1)
    expect_s3_class(plotNullDistribution(nd), "ggplot")
})

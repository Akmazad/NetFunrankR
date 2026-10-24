data(COPD_proteinCodingGenes, package = "NetFunrankR")

test_that("gene parsing handles pasted text and files", {
    expect_identical(.parseGenes("TP53, MDM2;EGFR\n\"KRAS\"  TP53"),
        c("TP53", "MDM2", "EGFR", "KRAS"))
    expect_identical(.parseGenes(NULL), character())
    f <- tempfile(fileext = ".csv")
    writeLines(c("TP53,1.2", "EGFR,0.4"), f)
    expect_true(all(c("TP53", "EGFR") %in% .readGeneFile(f)))
})

test_that("generated R code is valid R", {
    cfg <- .appConfig(NULL, NULL, NULL, "COPD", "Genes", 1L, 100)
    spec <- list(network = names(cfg$networks)[1],
        library = names(cfg$geneSets)[1], minSize = 5, maxSize = 500,
        alpha = 0.05, beyond = "linker")
    code <- .rCode(spec, c("TP53", "EGFR"), cfg)
    expect_silent(parse(text = code))
    expect_true(any(grepl("stringNetwork", code)))
    spec$network <- .upload
    expect_true(any(grepl("buildNetwork", .rCode(spec, "TP53", cfg))))
})

test_that("background computation returns a result", {
    .setAppConfig(.appConfig(NULL, NULL, NULL, "COPD", "Genes", 1L, 100))
    cfg <- .getAppConfig()
    spec <- list(genes = COPD_proteinCodingGenes[1:60],
        network = names(cfg$networks)[1], library = names(cfg$geneSets)[1],
        minSize = 5, maxSize = 500, alpha = 0.05, beyond = "linker",
        background = NULL)
    out <- .appCompute(spec)
    expect_s4_class(out$result, "NFResult")
    expect_true(out$seconds >= 0)
})

test_that("the dashboard accepts custom libraries and networks", {
    net <- toyNetwork()
    set.seed(11)
    gs <- list(A = sample(nodeNames(net), 30), B = sample(nodeNames(net), 30))
    cfg <- .appConfig(geneSets = list("Toy sets" = gs),
        networks = list("Toy network" = function() net),
        exampleGenes = character(0), exampleLabel = "none",
        geneLabel = "Toy IDs", workers = 1L, maxUploadMB = 10)
    expect_false(cfg$defaults)
    expect_false(.hasExample(cfg))
    .setAppConfig(cfg)
    out <- .appCompute(list(genes = c(gs$A[1:5], sample(nodeNames(net), 10)),
        network = "Toy network", library = "Toy sets", minSize = 1,
        maxSize = 500, alpha = 0.05, beyond = "linker", background = NULL))
    expect_setequal(enrichResults(out$result)$set, c("A", "B"))
    expect_true(any(grepl("your network", .rCode(out$spec, "G1", cfg))))
    expect_error(.appConfig(list(1), NULL, NULL, "x", "y", 1L, 1), "names")
    .setAppConfig(.appConfig(NULL, NULL, NULL, "COPD", "Genes", 1L, 100))
})

test_that("dashboard UI builds and the server validates inputs", {
    for (p in .appPackages) skip_if_not_installed(p)
    cfg <- .appConfig(NULL, NULL, NULL, "COPD", "Gene symbols", 1L, 100)
    .setAppConfig(cfg)
    expect_s3_class(.appUI(cfg), "shiny.tag.list")
    expect_s3_class(netFunrankRApp(), "shiny.appobj")
    shiny::testServer(.appServer, {
        session$setInputs(genes = "CHRNA3 CHRNA5 NOT_A_GENE",
            network = names(cfg$networks)[1],
            library = names(cfg$geneSets)[1], alpha = 0.05)
        expect_identical(st$genes(), c("CHRNA3", "CHRNA5", "NOT_A_GENE"))
        expect_true(st$iv$is_valid())
        session$setInputs(genes = "TP53")
        expect_false(st$iv$is_valid())
    })
})

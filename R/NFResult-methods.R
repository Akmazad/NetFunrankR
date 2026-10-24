#' @rdname NFResult-accessors
setMethod("enrichResults", "NFResult", function(x, ...) x@results)

#' @rdname NFResult-accessors
setMethod("queryGenes", "NFResult", function(x, ...) x@query)

#' @rdname NFResult-accessors
setMethod("unmatchedGenes", "NFResult", function(x, ...) x@unmatched)

#' @rdname NFResult-accessors
setMethod("nfParams", "NFResult", function(x, ...) x@params)

#' @rdname NFResult-accessors
setMethod("testedGeneSets", "NFResult", function(x, ...) x@geneSets)

#' @rdname evidenceClass
setMethod("evidenceClass", "NFResult", function(x, alpha = 0.05, ...) {
    r <- x@results
    o <- r$p_overlap < alpha
    b <- if (x@params$beyond == "distance") r$p_beyond < alpha else
        r$p_linker < alpha
    out <- ifelse(o & b, "overlap + network",
        ifelse(o, "overlap only", ifelse(b, "network only", "joint")))
    out[r$fdr_combined >= alpha] <- "none"
    factor(out, levels = .evidenceLevels)
})

.evidenceLevels <- c("overlap + network", "overlap only", "network only",
    "joint", "none")

#' @rdname NFResult-accessors
#' @param object an \code{\linkS4class{NFResult}}.
setMethod("show", "NFResult", function(object) {
    r <- object@results
    ev <- table(evidenceClass(object))
    cat("NFResult: decomposed network-proximity enrichment\n")
    cat(sprintf("  query: %d genes in network, %d not found\n",
        length(object@query), length(object@unmatched)))
    cat(sprintf("  network: %s\n", object@params$network))
    cat(sprintf("  %d gene sets tested; %d significant (FDR < 0.05)\n",
        nrow(r), sum(r$fdr_combined < 0.05)))
    cat(sprintf("    overlap + network %d | overlap only %d | ",
        ev[[1]], ev[[2]]))
    cat(sprintf("network only %d | joint %d\n", ev[[3]], ev[[4]]))
    top <- utils::head(r[r$fdr_combined < 0.05, c("set", "overlap",
        "n_linker", "fdr_combined")], 5)
    if (nrow(top)) {
        cat("  top gene sets:\n")
        top$fdr_combined <- signif(top$fdr_combined, 2)
        print(as.data.frame(top), row.names = FALSE)
    }
    invisible(object)
})

#' Coerce an NFResult to a data.frame
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param row.names,optional ignored.
#' @param alpha threshold used for the added \code{evidence} column.
#' @param ... ignored.
#' @returns \code{data.frame} of the results with an \code{evidence} column.
#' @exportS3Method base::as.data.frame
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' head(as.data.frame(res)[, c("set", "p_combined", "evidence")])
as.data.frame.NFResult <- function(x, row.names = NULL, optional = FALSE,
    alpha = 0.05, ...) {
    out <- as.data.frame(x@results)
    out$evidence <- as.character(evidenceClass(x, alpha))
    out
}

#' @rdname as.data.frame.NFResult
#' @exportMethod as.data.frame
setMethod("as.data.frame", "NFResult", as.data.frame.NFResult)

#' Write NetFunrankR results to a file
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param file output path; \code{.csv} or \code{.tsv}.
#' @param alpha threshold for the \code{evidence} column.
#' @returns \code{file}, invisibly.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' f <- exportResults(res, tempfile(fileext = ".csv"))
#' file.exists(f)
exportResults <- function(x, file, alpha = 0.05) {
    df <- as.data.frame(x, alpha = alpha)
    sep <- if (grepl("\\.tsv$", file, ignore.case = TRUE)) "\t" else ","
    utils::write.table(df, file, sep = sep, row.names = FALSE, quote = TRUE)
    invisible(file)
}

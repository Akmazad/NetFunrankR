#' Load a built-in pathway library as a named gene-set list
#'
#' @param db one of "KEGG", "Reactome", "WikiPathways".
#' @returns named list of gene symbol vectors.
#' @export
load_pathways <- function(db = c("KEGG", "Reactome", "WikiPathways")) {
  db <- match.arg(db)
  obj <- switch(db,
                KEGG = "KEGG_2021_Human",
                Reactome = "Reactome_Pathways_2024",
                WikiPathways = "WikiPathways_2024_Human")
  env <- new.env()
  utils::data(list = obj, package = "NetFunrankR", envir = env)
  as_geneset_list(get(obj, envir = env))
}

#' Built-in STRING v12 interactome
#'
#' @param evidence "noKB" (default) recomputes the STRING combined score
#'   without the curated-database and text-mining channels, so that the
#'   network is not derived from the pathway annotations being tested;
#'   "all" uses the published combined score.
#' @param score_th minimum score (0-1000). Defaults: 400 for "noKB",
#'   700 for "all" (similar network size, ~16,000 proteins).
#' @param ... passed to \code{\link{nf_network}}.
#' @returns an \code{nf_network}.
#' @export
string_network <- function(evidence = c("noKB", "all"), score_th = NULL, ...) {
  evidence <- match.arg(evidence)
  env <- new.env()
  utils::data("STRING_v12_PPI", package = "NetFunrankR", envir = env)
  ppi <- env$STRING_v12_PPI
  col <- if (evidence == "noKB") "nokb_score" else "combined_score"
  if (is.null(score_th)) score_th <- if (evidence == "noKB") 400 else 700
  if (score_th < 400 && evidence == "noKB")
    warning("bundled noKB edges are only complete for score >= 400")
  nf_network(ppi[ppi[[col]] >= score_th, c("from", "to")], ...)
}

#' Over-representation analysis (hypergeometric test)
#'
#' One-sided test \eqn{P(X \ge x)} for the overlap between the query genes
#' and every gene set, against a common background.
#'
#' @param query character vector of query genes.
#' @param genesets named list of gene vectors or an annotation data.frame.
#' @param universe background genes; default is every gene annotated in
#'   \code{genesets} (as in clusterProfiler/DAVID).
#' @param min_size,max_size size limits (genes in the universe).
#' @param p_adjust multiple-testing method.
#' @returns data.frame ordered by p-value.
#' @export
nf_ora <- function(query, genesets, universe = NULL,
                   min_size = 5L, max_size = Inf, p_adjust = "BH") {
  genesets <- as_geneset_list(genesets)
  if (is.null(universe)) universe <- unique(unlist(genesets, use.names = FALSE))
  q <- intersect(unique(query), universe)
  N <- length(universe)
  n <- length(q)
  res <- lapply(names(genesets), function(s) {
    Y <- intersect(genesets[[s]], universe)
    M <- length(Y)
    ov <- intersect(q, Y)
    x <- length(ov)
    data.frame(set = s, set_size = M, overlap = x,
               expected = n * M / N,
               gene_ratio = if (M > 0) x / M else NA_real_,
               pvalue = stats::phyper(x - 1L, M, N - M, n, lower.tail = FALSE),
               overlap_genes = paste(ov, collapse = "/"),
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res <- res[res$set_size >= min_size & res$set_size <= max_size, ]
  res$p.adjust <- stats::p.adjust(res$pvalue, method = p_adjust)
  res <- res[order(res$pvalue), ]
  rownames(res) <- NULL
  res
}

#' Exact null distribution of the proximity statistics for one gene set
#'
#' Returns the probability mass functions used by
#' \code{\link{nf_proximity_test}} for a single set, e.g. for plotting.
#'
#' @inheritParams nf_proximity_test
#' @param set name of the gene set.
#' @returns list with observed values and pmfs of the overlap (\code{O}),
#'   total distance (\code{S}) and beyond-overlap distance (\code{Sb}).
#' @export
nf_null_distribution <- function(query, genesets, net, set) {
  genesets <- as_geneset_list(genesets)[set]
  d <- set_distances(net, genesets)[, 1]
  qi <- which(net$genes %in% query)
  B <- max(net$bin)
  cnt <- matrix(tabulate(d * B + net$bin, nbins = B * (max(d) + 1L)), nrow = B)
  f <- cnt / pmax(rowSums(cnt), 1)
  member <- d[qi] == 0L
  pm <- exact_null_cpp(f, tabulate(net$bin[qi], nbins = B),
                       tabulate(net$bin[qi][!member], nbins = B))
  c(pm, list(n = length(qi), overlap = sum(member), S_obs = sum(d[qi]),
             d_query = stats::setNames(d[qi], net$genes[qi])))
}

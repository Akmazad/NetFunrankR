#' Build a NetFunrankR interactome object
#'
#' Converts an undirected edge list into the network object used by all
#' NetFunrankR tests: the largest connected component (LCC) of the simple
#' graph, its sparse adjacency matrix, node degrees and degree bins.
#'
#' Degree bins follow Guney et al. (2016): nodes are sorted by degree and
#' consecutive degree values are pooled until a bin holds at least
#' \code{min_bin_size} nodes. Random replacements for a query gene are drawn
#' from its own bin, so the null preserves the degree profile of the query.
#'
#' @param edges data.frame/matrix whose first two columns are node IDs.
#' @param lcc logical; restrict to the largest connected component
#'   (default \code{TRUE}; required for finite shortest-path distances).
#' @param min_bin_size minimum number of nodes per degree bin (default 100).
#'
#' @returns An object of class \code{nf_network}.
#' @export
nf_network <- function(edges, lcc = TRUE, min_bin_size = 100L) {
  edges <- as.data.frame(edges)[, 1:2]
  colnames(edges) <- c("from", "to")
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges <- edges[!is.na(edges$from) & !is.na(edges$to) &
                   edges$from != edges$to, , drop = FALSE]

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  g <- igraph::simplify(g)
  if (lcc) {
    comp <- igraph::components(g)
    g <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  }

  genes <- igraph::V(g)$name
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  A <- methods::as(A, "CsparseMatrix")
  deg <- as.integer(Matrix::rowSums(A))
  names(deg) <- genes

  structure(
    list(graph = g,
         genes = genes,
         adj = A,
         degree = deg,
         bin = degree_bins(deg, min_bin_size),
         min_bin_size = min_bin_size),
    class = "nf_network"
  )
}

#' @export
print.nf_network <- function(x, ...) {
  cat(sprintf("NetFunrankR network: %d nodes, %d edges, %d degree bins (min size %d)\n",
              length(x$genes), igraph::ecount(x$graph),
              length(unique(x$bin)), x$min_bin_size))
  invisible(x)
}

#' Degree binning (Guney et al. 2016)
#'
#' @param deg integer vector of node degrees.
#' @param min_bin_size minimum number of nodes per bin.
#' @returns integer vector of bin labels (1 = lowest degree), same order as
#'   \code{deg}.
#' @export
degree_bins <- function(deg, min_bin_size = 100L) {
  tab <- table(deg)
  dvals <- as.integer(names(tab))
  counts <- as.integer(tab)
  bin_of_value <- integer(length(dvals))
  b <- 1L
  acc <- 0L
  for (i in seq_along(dvals)) {
    bin_of_value[i] <- b
    acc <- acc + counts[i]
    if (acc >= min_bin_size && i < length(dvals)) {
      b <- b + 1L
      acc <- 0L
    }
  }
  # merge an undersized last bin into its predecessor
  last <- bin_of_value == b
  if (b > 1L && sum(counts[last]) < min_bin_size)
    bin_of_value[last] <- b - 1L
  bin_of_value[match(deg, dvals)]
}

#' Convert gene sets to a named list restricted to valid IDs
#'
#' Accepts a named list, or a two-column data.frame such as the bundled
#' Enrichr libraries (\code{GeneSymbol}, \code{PathwayName}).
#'
#' @param genesets named list of character vectors or a data.frame.
#' @returns named list of unique character vectors.
#' @export
as_geneset_list <- function(genesets) {
  if (is.data.frame(genesets)) {
    gcol <- intersect(c("GeneSymbol", "gene", "Gene"), colnames(genesets))[1]
    scol <- intersect(c("PathwayName", "term", "Term", "set"), colnames(genesets))[1]
    if (is.na(gcol) || is.na(scol))
      stop("data.frame gene sets need GeneSymbol/PathwayName columns")
    genesets <- split(as.character(genesets[[gcol]]),
                      as.character(genesets[[scol]]))
  }
  if (!is.list(genesets) || is.null(names(genesets)))
    stop("genesets must be a named list or an annotation data.frame")
  lapply(genesets, function(x) unique(as.character(x[!is.na(x)])))
}

#' Shortest-path distance from every node to every gene set
#'
#' Computes \eqn{d_Y(v) = \min_{y \in Y} d(v, y)} for all nodes \eqn{v} and
#' all gene sets \eqn{Y} simultaneously, by a multi-source breadth-first
#' search expressed as repeated sparse matrix products (one product per BFS
#' layer, shared by all sets).
#'
#' @param net an \code{nf_network}.
#' @param genesets named list of gene vectors.
#' @returns integer matrix (nodes x sets); \code{NA} for sets without any
#'   gene in the network.
#' @export
set_distances <- function(net, genesets) {
  genes <- net$genes
  idx <- lapply(genesets, function(s) which(genes %in% s))
  K <- length(idx)
  V <- length(genes)
  D <- matrix(NA_integer_, V, K, dimnames = list(genes, names(genesets)))
  frontier <- matrix(FALSE, V, K)
  for (k in seq_len(K)) frontier[idx[[k]], k] <- TRUE
  D[frontier] <- 0L
  visited <- frontier
  step <- 0L
  while (any(frontier)) {
    step <- step + 1L
    nb <- as.matrix(net$adj %*% (frontier * 1)) > 0
    frontier <- nb & !visited
    D[frontier] <- step
    visited <- visited | frontier
  }
  D
}

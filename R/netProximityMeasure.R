#' Network proximity between two gene modules (NetFunrankR v1 interface)
#'
#' Symmetric module-to-module proximity
#' \eqn{d_{XY} = (\sum_x \min_y d(x,y) + \sum_y \min_x d(x,y)) / (|X| + |Y|)}
#' with a degree-binned permutation null for the query module. Kept for
#' backward compatibility; \code{\link{nf_proximity_test}} is recommended.
#'
#' @param net an \code{nf_network}, or an igraph object.
#' @param module_a_genes query (disease) genes.
#' @param module_b_genes gene set (pathway) genes.
#' @param nPermute number of permutations.
#' @returns numeric vector: observed distance, z-score, empirical p-value.
#' @export
network.proximity <- function(net, module_a_genes, module_b_genes, nPermute = 1000) {
  if (inherits(net, "igraph")) net <- nf_network(igraph::as_data_frame(net)[, 1:2])
  r <- nf_proximity_perm(module_a_genes, list(Y = module_b_genes), net,
                         n_perm = nPermute, statistic = "symmetric", min_size = 1L)
  c(observed = r$observed, z = r$z, p = r$p_perm)
}

#' Network proximity score on the built-in STRING interactome (v1 interface)
#'
#' @param string_PPI_score_th STRING score threshold (noKB evidence).
#' @param module_a query (disease) genes.
#' @param module_b gene set (pathway) genes.
#' @param nPermute number of permutations.
#' @returns see \code{\link{network.proximity}}.
#' @export
quantify_network_proximity_score <- function(string_PPI_score_th = 400,
                                             module_a = c(), module_b = c(),
                                             nPermute = 1000) {
  network.proximity(string_network("noKB", string_PPI_score_th),
                    module_a, module_b, nPermute)
}

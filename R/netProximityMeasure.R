#' Title Compute network proximity scores between disease genes and the genes in a pathway underlying PPI network (after filtering the PPI)
#'
#' @param string_PPI_score_th a threshold for filtering StringDB PPI
#' @param module_a a vector of disease genes
#' @param module_b a vector of pathway genes
#' @param nPermute number of permutation for significance calculation
#'
#' @returns returns the proximity scores
#' @export
#'
#' @examples
quantify_network_proximity_score <- function(
    string_PPI_score_th,
    module_a = c(),
    module_b = c(),
    nPermute
){

  data(protein.links)
  # extract PPI for Human proteins that are enriched in two pathways
  string.ppi.df = protein.links %>%
    dplyr::filter(combined_score > string_PPI_score_th) %>%
    dplyr::select(-combined_score) %>%
    # CoDiNA::OrderNames() %>%
    # unique() %>%
    igraph::graph_from_data_frame(directed = F)

  # browser()
  # get network proximity score
  network.proximity(net = string.ppi.df,
                    module_a_genes = module_a,
                    module_b_genes = module_b,
                    nPermute = nPermute)
}

#' Title Compute network proximity scores between disease genes and the genes in a pathway underlying a PPI network
#'
#' @param net an igraph object for the PPI network
#' @param module_a_genes a vector of disease genes
#' @param module_b_genes a vector of pathway genes
#' @param nPermute number of permutation for significance calculation
#'
#' @returns returns the proximity scores
#' @export
#'
#' @examples
network.proximity <- function(net,
                              module_a_genes,
                              module_b_genes,
                              nPermute){


  d <- module_a_genes
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  # d_temp <- unique(d[keep])
  # # get PPI partners of to expand the space
  # d <- c(d_temp,
  #        getPPIpartners(ppi.net = net,
  #                       geneList = d_temp,
  #                       hop = neighbourhood_th)) %>% unique()

  t <- module_b_genes
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])
  # t_temp <- unique(t[keep])
  # # get PPI partners of to expand the space
  # t <- c(t_temp,
  #        getPPIpartners(ppi.net = net,
  #                       geneList = t_temp,
  #                       hop = neighbourhood_th)) %>% unique()

  # browser()
  d_td <- get_d(shortest.paths(net, v = t, to=d))
  # return(paste0("[Observed distance: ", d_td))

  z <- permuteTest(net, t, d, d_td, nPermute)
  # One-sided (left-tail) hypothesis testing [alternate hypo: observed z-score < mean of random z-scores]
  p <- pnorm(z)

  # return(paste0("[Observed distance: ", d_td,"; z-score: ", z, "; p-value: ", p,"]"))
  return(c(d_td, z, p))
}

get_d <- function(dist_matrix){
    clean_matrix <- dist_matrix %>%
    as.data.frame() %>%
    filter(rowSums(is.infinite(as.matrix(.))) != ncol(.)) %>%
    dplyr::select(where(~ !all(is.infinite(.x)))) %>%
    as.matrix()

  d <- mean(c(apply(clean_matrix,1,min), apply(clean_matrix,2,min)))

  return(d)
}
getRandD <- function(a_degree, net){
  return(base::sample(which(igraph::degree(net) == a_degree),1, replace = F))
}
permuteTest <- function(net, t, d, d_td, N){
  r <- c()
  for (i in 1:N) {
    t_rand <- sapply(as.numeric(igraph::degree(net, v = t)), getRandD,net)
    d_rand <- sapply(as.numeric(igraph::degree(net, v = d)), getRandD,net)
    d_td_rand <- get_d(shortest.paths(net, v = t_rand %>% unique(), to=d_rand %>% unique()))
    r<-c(r,d_td_rand)
  }
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m)/s
  return(z)
}

#' Title
#'
#' @param neighbourhood_th
#' @param string_PPI_score_th
#' @param module_a
#' @param module_b
#' @param nPermute
#'
#' @returns
#' @export
#'
#' @examples
quantify_network_proximity_score <- function(
    neighbourhood_th,
    string_PPI_score_th,
    module_a = c(),
    module_b = c(),
    nPermute
){

  data(protein.links)
  # extract PPI for Human proteins that are enriched in two pathways
  string.ppi.df = protein.links %>%
    dplyr::filter(combined_score > string_PPI_score_th) %>%
    # dplyr::select(-combined_score)
    CoDiNA::OrderNames() %>%
    unique() %>%
    igraph::graph_from_data_frame(directed = F)

  # get network proximity score
  network.proximity(net = string.ppi.df, module_a_genes = module_a,
                    module_b_genes = module_b, nPermute = nPermute)
}

#' Title
#'
#' @param net
#' @param module_a_genes
#' @param module_b_genes
#' @param nPermute
#'
#' @returns
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

  t <- module_b_genes
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])

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
    select(where(~ !all(is.infinite(.x)))) %>%
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

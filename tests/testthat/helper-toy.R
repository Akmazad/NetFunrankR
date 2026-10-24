toyNetwork <- function(n = 300, m = 3, seed = 1) {
    set.seed(seed)
    g <- igraph::sample_pa(n, m = m, directed = FALSE)
    igraph::V(g)$name <- paste0("G", seq_len(n))
    buildNetwork(g, minBinSize = 20L)
}

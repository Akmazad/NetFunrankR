suppressPackageStartupMessages({
  library(igraph)
  library(SummarizedExperiment)
})

# ------------------------------------------------------------
# 1) Build a weighted probabilistic PPI graph
#    Assumes STRING-style table with columns: from, to, combined_score
#    IMPORTANT: from/to must already be in the same ID space as rownames(se)
#               and the genes in gs (preferably Entrez IDs).
# ------------------------------------------------------------
ppi_df <- protein.links
build_probabilistic_ppi_weighted <- function(ppi_df,
                                             score_col = "combined_score",
                                             from_col  = "from",
                                             to_col    = "to",
                                             min_score = 400) {
  stopifnot(all(c(from_col, to_col, score_col) %in% colnames(ppi_df)))

  df <- ppi_df[!is.na(ppi_df[[score_col]]) & ppi_df[[score_col]] >= min_score, ,
               drop = FALSE]

  df[[from_col]] <- as.character(df[[from_col]])
  df[[to_col]]   <- as.character(df[[to_col]])
  df[[score_col]] <- as.numeric(df[[score_col]])

  # STRING confidence -> probability in (0,1)
  max_sc <- max(df[[score_col]], na.rm = TRUE)
  prob   <- df[[score_col]] / max_sc
  prob   <- pmin(pmax(prob, 1e-6), 0.999999)

  # Stable positive edge cost: high confidence => low cost
  weight <- (1 - prob) + 1e-6

  g <- graph_from_data_frame(
    data.frame(
      from   = df[[from_col]],
      to     = df[[to_col]],
      prob   = prob,
      weight = weight,
      stringsAsFactors = FALSE
    ),
    directed = FALSE
  )

  E(g)$prob   <- weight <- E(g)$prob %||% prob
  E(g)$weight <- weight

  g
}

# small infix helper
`%||%` <- function(x, y) if (is.null(x)) y else x


# ------------------------------------------------------------
# 2) Multi-source shortest path from one source set to all nodes
# ------------------------------------------------------------
compute_multi_source_dist <- function(g, sources, large_dist = 10) {
  nodes   <- V(g)$name
  sources <- intersect(as.character(sources), nodes)

  if (length(sources) == 0L) {
    out <- rep(large_dist, length(nodes))
    names(out) <- nodes
    return(out)
  }

  dmat <- distances(
    g,
    v       = sources,
    to      = nodes,
    mode    = "all",
    weights = E(g)$weight
  )

  dmin <- apply(dmat, 2, min)
  dmin[is.infinite(dmin)] <- large_dist
  names(dmin) <- nodes
  dmin
}


# ------------------------------------------------------------
# 3) Original symmetric proximity:
#    d_{X,Y} = 1/(|X|+|Y|) [ sum_x min_y d(x,y) + sum_y min_x d(x,y) ]
#    implemented via precomputed d_X(.) and d_Y(.)
# ------------------------------------------------------------
compute_symmetric_proximity <- function(X, Y, d_X, d_Y, large_dist = 10) {
  X <- intersect(as.character(X), names(d_Y))
  Y <- intersect(as.character(Y), names(d_X))

  if (length(X) == 0L || length(Y) == 0L) {
    return(large_dist)
  }

  d_X_to_Y <- mean(d_Y[X], na.rm = TRUE)  # mean_x min_y d(x,y)
  d_Y_to_X <- mean(d_X[Y], na.rm = TRUE)  # mean_y min_x d(x,y)

  ((length(X) * d_X_to_Y) + (length(Y) * d_Y_to_X)) / (length(X) + length(Y))
}


# ------------------------------------------------------------
# 4) Query gene extractor from se
#    This is the one part you may want to swap depending on your method:
#    - "de": adjusted p-value threshold
#    - "top_n": top |STAT| genes
# ------------------------------------------------------------
extract_query_genes <- function(se,
                                mode = c("de", "top_n"),
                                alpha = 0.05,
                                top_n = 200,
                                stat_col = "limma.STAT",
                                padj_col = "ADJ.PVAL") {
  mode <- match.arg(mode)

  rr <- as.data.frame(rowData(se))
  gene_ids <- rownames(se)

  if (is.null(gene_ids) || length(gene_ids) == 0L) {
    stop("rownames(se) must contain gene IDs.")
  }

  if (mode == "de") {
    if (!padj_col %in% colnames(rr)) {
      stop(sprintf("rowData(se) must contain '%s' for mode='de'.", padj_col))
    }
    keep <- !is.na(rr[[padj_col]]) & rr[[padj_col]] <= alpha
    q <- gene_ids[keep]
  } else {
    if (!stat_col %in% colnames(rr)) {
      stop(sprintf("rowData(se) must contain '%s' for mode='top_n'.", stat_col))
    }
    ord <- order(abs(rr[[stat_col]]), decreasing = TRUE, na.last = NA)
    q <- gene_ids[head(ord, top_n)]
  }

  unique(as.character(q))
}


# ------------------------------------------------------------
# 5) Factory: returns a benchmarkable method(se, gs)
#    This is the function to pass to runEA()
# ------------------------------------------------------------
make_netfunrankr_benchmark <- function(ppi_df,
                                       score_col = "combined_score",
                                       from_col  = "from",
                                       to_col    = "to",
                                       min_score = 400,
                                       n_perm    = 200,
                                       large_dist = 10,
                                       query_mode = c("de", "top_n"),
                                       alpha = 0.05,
                                       top_n = 200,
                                       stat_col = "limma.STAT",
                                       padj_col = "ADJ.PVAL",
                                       restrict_to_gs = TRUE,
                                       seed = 1L) {

  query_mode <- match.arg(query_mode)

  # Build graph once
  g_full <- build_probabilistic_ppi_weighted(
    ppi_df    = ppi_df,
    score_col = score_col,
    from_col  = from_col,
    to_col    = to_col,
    min_score = min_score
  )

  function(se, gs) {
    set.seed(seed)

    # gs is expected to be a named list: pathway_id -> gene vector
    if (!is.list(gs) || is.null(names(gs))) {
      stop("gs must be a named list of gene sets.")
    }

    # Harmonize gene sets
    gs <- lapply(gs, function(x) unique(as.character(x)))
    gs <- gs[lengths(gs) > 0L]

    # Extract query genes from dataset
    X <- extract_query_genes(
      se       = se,
      mode     = query_mode,
      alpha    = alpha,
      top_n    = top_n,
      stat_col = stat_col,
      padj_col = padj_col
    )

    if (length(X) == 0L) {
      # conservative fallback: all 1s
      p <- rep(1, length(gs))
      names(p) <- names(gs)
      return(p)
    }

    # Restrict graph for speed
    g <- g_full
    if (restrict_to_gs) {
      needed <- unique(c(X, unlist(gs, use.names = FALSE)))
      needed <- intersect(needed, V(g_full)$name)
      g <- induced_subgraph(g_full, vids = needed)
    }

    nodes <- V(g)$name
    X     <- intersect(X, nodes)

    if (length(X) == 0L) {
      p <- rep(1, length(gs))
      names(p) <- names(gs)
      return(p)
    }

    # Keep only gene sets that intersect graph
    gs_in <- lapply(gs, intersect, y = nodes)
    names(gs_in) <- names(gs)

    # Precompute d_X(.) once
    d_X <- compute_multi_source_dist(g, X, large_dist = large_dist)

    # Precompute d_Y(.) once per pathway and observed score
    path_names <- names(gs_in)
    K <- length(path_names)

    dY_list  <- vector("list", K)
    obs_dist <- rep(large_dist, K)

    for (k in seq_len(K)) {
      Y <- gs_in[[k]]
      if (length(Y) == 0L) {
        dY_list[[k]] <- NULL
        obs_dist[k]  <- large_dist
        next
      }

      d_Y <- compute_multi_source_dist(g, Y, large_dist = large_dist)
      dY_list[[k]] <- d_Y
      obs_dist[k]  <- compute_symmetric_proximity(
        X = X, Y = Y, d_X = d_X, d_Y = d_Y, large_dist = large_dist
      )
    }

    # Permutation null
    null_mat <- matrix(NA_real_, nrow = n_perm, ncol = K)
    colnames(null_mat) <- path_names

    for (b in seq_len(n_perm)) {
      Rb   <- sample(nodes, length(X), replace = FALSE)
      d_Rb <- compute_multi_source_dist(g, Rb, large_dist = large_dist)

      for (k in seq_len(K)) {
        Y <- gs_in[[k]]
        if (length(Y) == 0L || is.null(dY_list[[k]])) {
          null_mat[b, k] <- large_dist
          next
        }

        null_mat[b, k] <- compute_symmetric_proximity(
          X = Rb,
          Y = Y,
          d_X = d_Rb,
          d_Y = dY_list[[k]],
          large_dist = large_dist
        )
      }
    }

    # Smaller distance = more enriched
    # permutation p-value with +1 correction (never zero)
    pvals <- vapply(seq_len(K), function(k) {
      obs  <- obs_dist[k]
      null <- null_mat[, k]
      (sum(null <= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
    }, numeric(1))

    names(pvals) <- path_names
    pvals
  }
}


library(GSEABenchmarkeR)
library(EnrichmentBrowser)

# 1) Load benchmark compendium
geo2kegg <- loadEData("geo2kegg")
geo2kegg <- maPreproc(geo2kegg[1:5])
geo2kegg <- runDE(geo2kegg, de.method = "limma", padj.method = "flexible")

# 2) KEGG gene sets
kegg.gs <- getGenesets(org = "hsa", db = "kegg")

# 3) Your PPI table: MUST be in Entrez IDs to match rownames(se) and kegg.gs
# Example expected columns: from, to, combined_score
# ppi_df <- readRDS("string_entrez_ppi.rds")

# 4) Create benchmark-ready method
netfun_prob <- make_netfunrankr_benchmark(
  ppi_df      = ppi_df,
  score_col   = "combined_score",
  from_col    = "from",
  to_col      = "to",
  min_score   = 700,        # 700 is often a good benchmark compromise
  n_perm      = 100,        # start with 100 for benchmarking
  query_mode  = "de",       # or "top_n"
  alpha       = 0.05,
  top_n       = 200,
  stat_col    = "limma.STAT",
  padj_col    = "ADJ.PVAL",
  restrict_to_gs = TRUE,
  seed        = 1L
)

# 5) Run your method together with comparators
res.dir <- tempdir()

# built-in methods
res.std <- runEA(
  geo2kegg,
  methods   = c("ora", "gsea", "camera", "spia"),
  gs        = kegg.gs,
  perm      = 0, , grn = my_grn,
  save2file = TRUE,
  out.dir   = res.dir
)

# custom method
method <- netfun_prob
res.custom <- runEA(
  geo2kegg,
  method,
  kegg.gs,
  save2file = TRUE,
  out.dir   = res.dir
)

# read saved results
ea.rtimes <- readResults(
  res.dir,
  names(geo2kegg),
  methods = c("ora", "gsea", "camera", "method", "spia"),
  type    = "runtime"
)

ea.ranks <- readResults(
  res.dir,
  names(geo2kegg),
  methods = c("ora", "gsea", "camera", "method", "spia"),
  type    = "ranking"
)


# plot the comparison (runtime and percentage of significant sets)
bpPlot(ea.rtimes, what = "runtime")
sig.sets <- evalNrSigSets(ea.ranks, alpha = 0.05, padj = "BH")
bpPlot(sig.sets, what = "sig.sets")

# my_grn <- EnrichmentBrowser::compileGRN(
#   org = "hsa",
#   db = "kegg",
#   act.inh = TRUE,
#   map2entrez = TRUE,
#   keep.type = FALSE,
#   kegg.native = FALSE
# )

# method_list <- c("ora", "gsea", "camera", "netfun_prob")
#
# res.all <- runEA(
#   geo2kegg,
#   methods   = method_list,
#   gs        = kegg.gs,
#   perm      = 0,
#   save2file = TRUE, grn = my_grn,
#   out.dir   = res.dir
# )

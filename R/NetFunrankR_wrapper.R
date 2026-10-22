#' Run the complete NetFunrankR analysis
#'
#' Runs the decomposed network-proximity test (\code{\link{nf_proximity_test}})
#' together with classical over-representation analysis, writes the result
#' table and (optionally) the figures to \code{outdir}.
#'
#' @param disease_genes query genes (HGNC symbols for the built-in data).
#' @param PE_ORA_pop pathway library: "KEGG", "Reactome", "WikiPathways", or a
#'   named list / annotation data.frame of gene sets.
#' @param network an \code{nf_network}; default \code{string_network("noKB")}.
#' @param min_size,max_size gene-set size limits (genes in the network).
#' @param alpha FDR threshold for plots and the evidence classes.
#' @param doPlot logical; write figures.
#' @param plot_width,plot_height figure size (inches).
#' @param outdir output directory.
#'
#' @returns (invisibly) the \code{nf_result} data.frame.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- NetFunrankR_wrapper(disease_genes = COPD_proteinCodingGenes)
#' head(res)
#' }
NetFunrankR_wrapper <- function(disease_genes = NULL,
                                PE_ORA_pop = "KEGG",
                                network = NULL,
                                min_size = 5L,
                                max_size = 500L,
                                alpha = 0.05,
                                doPlot = TRUE,
                                plot_width = 10,
                                plot_height = 8,
                                outdir = "./NetFunrankR_results/") {
  if (is.null(disease_genes) || length(disease_genes) == 0)
    stop("Disease genes are required for this analysis")
  gs <- if (is.character(PE_ORA_pop) && length(PE_ORA_pop) == 1L)
    load_pathways(PE_ORA_pop) else as_geneset_list(PE_ORA_pop)
  lib <- if (is.character(PE_ORA_pop) && length(PE_ORA_pop) == 1L) PE_ORA_pop else "custom"
  if (is.null(network)) network <- string_network("noKB")

  res <- nf_proximity_test(disease_genes, gs, network,
                           min_size = min_size, max_size = max_size)
  res$evidence <- nf_evidence_class(res, alpha)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- as.data.frame(res)
  data.table::fwrite(out, file.path(outdir, paste0(lib, "_NetFunrankR_result.csv")))
  message(sprintf("%d query genes in network (%d not found); %d sets tested; %d sets with FDR(combined) < %.2f",
                  res$n_query_net[1], length(attr(res, "query_not_in_network")),
                  nrow(res), sum(res$fdr_combined < alpha), alpha))

  if (doPlot) {
    save_plot <- function(p, name) {
      ggplot2::ggsave(file.path(outdir, paste0(lib, "_", name, ".pdf")), p,
                      width = plot_width, height = plot_height, limitsize = FALSE)
      saveRDS(p, file.path(outdir, paste0(lib, "_", name, ".rds")))
    }
    save_plot(plot_evidence_map(res, alpha), "evidence_map")

    sig_ora <- res$set[res$fdr_hyper < alpha]
    sig_net <- res$set[res$fdr_combined < alpha]
    save_plot(plot_venn(sig_ora, sig_net,
                        c("Over-representation", "NetFunrankR")), "venn")
    common <- res[res$set %in% intersect(sig_ora, sig_net), ]
    if (nrow(common) > 0)
      save_plot(plot_double_bar(data.frame(set = common$set, p1 = common$p_hyper,
                                           p2 = common$p_combined),
                                c("Over-representation", "NetFunrankR")), "common_sets")
    only_net <- res[res$set %in% setdiff(sig_net, sig_ora), ]
    if (nrow(only_net) > 0)
      save_plot(plot_lollipop(data.frame(set = only_net$set, pvalue = only_net$p_combined)),
                "NetFunrankR_only_sets")
    top_net <- res[res$evidence %in% c("network only", "overlap + network"), ]
    if (nrow(top_net) > 0) {
      s <- top_net$set[which.min(top_net$p_linker)]
      save_plot(plot_linker_network(disease_genes, gs[[s]], network) +
                  ggplot2::ggtitle(s), "linker_network")
    }
  }
  invisible(res)
}

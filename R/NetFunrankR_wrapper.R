#' Title This function works as a wrapper for running the whole functional analysis using network medicine based approach.
#'
#' @param disease_genes Input genes (a vector) that are used functional enrichment analysis. Default is "COPD-related genes".
#' @param PathwayRanker_Type A string that indicates the type of enrichment to be done. Options are: "ORA" (Hypergeometric-distribution based over-representation test), and "TW" (Topology-aware network medicine-based approach), "Both" (Default: both ORA and TW so that the results can be compared).
#' @param PE_ORA_pop A string that indicates the set of known annotaiton of pathways or functional term. Current version has the following options: "KEGG" (default), "Reactome" and "WikiPathway".
#' @param nPermute A number (default: 1000) to indicate the number of permutation iteration to be conducted for hypothesis testing in Network-medicine based functional enrichment
#' @param string_PPI_score_th A number (default: 900) to filter StringDB PPIs only above the given value.
#' @param doPlot A logical value (default: TRUE) to specify whether to plot the results.
#' @param plot_width A number (default: 10) to specify the output plot width.
#' @param plot_height A number (default: 10) to specify the output plot height.
#' @param outdir A string to indicate the output directory to save all the results.
#'
#' @returns ""
#' @export ""
#'
#' @examples
#' #' library(NetFunrankR)
#' library(dplyr)
#' NetFunrankR_wrapper(disease_genes = COPD_proteinCodingGenes,
#'               doPlot = TRUE,
#'               nPermute = 100)
NetFunrankR_wrapper <- function(
    disease_genes = NULL,
    PathwayRanker_Type = "Both", # "c('ORA','TW', "Both")
    PE_ORA_pop = "KEGG", # c('KEGG', 'Reactome', 'WikiPathway')
    # PE_ORA_cutoff_col = "p.adjust", # c("p.adjust", "pvalue")
    # PE_ORA_cutoff_th = 0.05,
    nPermute = 10, # number of permutation iteration for hypothesis testing
    string_PPI_score_th = 900, # controlling PPI confidence
    doPlot = TRUE,
    plot_width = 10,
    plot_height = 8,
    outdir = "./NetFunrankR_results/"
){
  if(is.null(disease_genes))
    stop("Disease genes are required for this analysis!!")
  if(nPermute < 10)
    stop("Permutation number is very small. Please consider increasing it (100 or more - larger, better, but it will run longer) !!")


  dir.create(outdir, recursive = T)

  # 1) source the ORA geneset -----------
  pathway.data = NULL
  if(PE_ORA_pop == "KEGG"){
    data(KEGG_2021_Human)
    pathway.data = KEGG_2021_Human
  }else if(PE_ORA_pop == "Reactome"){
    data(Reactome_Pathways_2024)
    pathway.data = Reactome_Pathways_2024
  }else{
    data(WikiPathways_2024_Human)
    pathway.data = WikiPathways_2024_Human
  }
  # ----

  # 2) Pathway Ranking -----
  rankedPathways = NULL
  if(PathwayRanker_Type == "ORA"){
    set.seed(123)

    ## 2.1) do Pathway Enrichment test (ORA) [output-file]
    rankedPathways <- PathwayEnrichment_ORA(givenSet = disease_genes,
                                            pathway.data = pathway.data,
                                            PE_ORA_pop = PE_ORA_pop)
    rankedPathways %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_ORA_result.csv"))
  }else if(PathwayRanker_Type == "TW"){
    set.seed(123)

    ## 2.2) Do Functional proximity [output-file]
    rankedPathways <- PathwayEnrichment_TopoAw(
      string_PPI_score_th = string_PPI_score_th,
      nPermute = nPermute,
      pathway.data = pathway.data,
      givenSet = disease_genes)

    rankedPathways %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_TopoAw_result.csv"))
  }else{
    set.seed(123)

    message(cat("\n===============================================================\n"))
    message(cat("[1] Over-representation (Hypergeometric dist.) is in progress: \n"))
    message(cat("===============================================================\n"))
    rankedPathways_ORA <- PathwayEnrichment_ORA(givenSet = disease_genes,
                                            pathway.data = pathway.data,
                                            PE_ORA_pop = PE_ORA_pop)
    rankedPathways_ORA %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_ORA_result.csv"))
    message(cat("\n======================      [DONE]     ========================\n"))

    message(cat("\n===============================================================\n"))
    message(cat("[2] Network-proxmity-based analysis is in progress: \n"))
    message(cat("=================================================================\n"))
    #
    set.seed(123)

    rankedPathways_TW <- PathwayEnrichment_TopoAw(
      string_PPI_score_th = string_PPI_score_th,
      nPermute = nPermute,
      pathway.data = pathway.data,
      givenSet = disease_genes)
    rankedPathways_TW %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_TopoAw_result.csv"))

    message(cat("\n======================      [DONE]     ========================\n"))

    # browser()
    # compare

    # remove later ----
    rankedPathways_ORA <- fread(paste0(outdir, PE_ORA_pop, "_PE_ORA_result.csv"))
    rankedPathways_TW <- fread(paste0(outdir, PE_ORA_pop, "_PE_TopoAw_result.csv"))
    rankedPathways_ORA_filt <- rankedPathways_ORA %>% dplyr::filter(pvalue < 0.05)
    rankedPathways_TW_filt <- rankedPathways_TW  %>% dplyr::filter(pvalue < 0.05)
    # ----
    if(doPlot){
      # plot common terms/pathways ---------
      # venn
      common.terms.vennplt.obj <- plot_venn(
        vec1 = rankedPathways_ORA_filt %>% dplyr::pull(Description),
        vec2 = rankedPathways_TW_filt %>% dplyr::pull(Description))
      save(common.terms.vennplt.obj,
           file = paste0(outdir, "CommonFunctions_ORA_TW.vennplt.obj.RData"))
      common.terms.vennplt.obj %>%
        ggsave(filename = paste0(outdir, "CommonFunctions_ORA_TW_vennplot.pdf"),
               width = plot_width, height = plot_height, limitsize = F)
      #
      plt.df <- dplyr::inner_join(
        rankedPathways_ORA_filt,
        rankedPathways_TW_filt,
        by = "Description")

      if(plt.df %>% nrow() > 1){
        common.terms.barplt.obj <- plt.df %>%
          plot_double_bar()
        # save the plot to the result directory
        save(common.terms.barplt.obj,
             file = paste0(outdir, "CommonFunctions_ORA_TW.barplt.obj.RData"))
        common.terms.barplt.obj %>%
          ggsave(filename = paste0(outdir, "CommonFunctions_ORA_TW_barplot.pdf"),
                 width = plot_width, height = plot_height, limitsize = F)
      }else
        warning(cat("No common terms found between ORA and TW-based analysis with given filter"))


      #-----

      # plot lollipop plots for distinct ones [left exclusive joins] ---------
      plt.df <- dplyr::anti_join(
        rankedPathways_ORA_filt,
        rankedPathways_TW_filt,
        by = c("Description" = "Description"))

      if(plt.df %>% nrow() > 1){
        only_ORA_terms.plt.obj <- plt.df %>%
          plot_lollipop(color_palette = c("#616530FF", "#dce58f"))
        save(only_ORA_terms.plt.obj,
             file = paste0(outdir, "Only_ORA_terms.plt.obj.RData"))
        only_ORA_terms.plt.obj %>%
          ggsave(filename = paste0(outdir, "only_ORA_terms_lollipop.pdf"),
                 width = plot_width, height = plot_height, limitsize = F)
      }else
        warning(cat("No exclusive terms found for ORA based analysis with given filter"))

      # plot lollipop plots for distinct ones [left exclusive joins] ---------
       plt.df <- dplyr::anti_join(
        rankedPathways_TW_filt,
        rankedPathways_ORA_filt,
        by = c("Description" = "Description"))

      if(plt.df %>% nrow() > 1){
        only_TW_terms.plt.obj <- plt.df %>%
          plot_lollipop(color_palette = c("#CC8214FF", "#f6e8b1"))
        save(only_TW_terms.plt.obj,
             file = paste0(outdir, "Only_TW_terms.plt.obj.RData"))
        only_TW_terms.plt.obj %>%
          ggsave(filename = paste0(outdir, "only_TW_terms_lollipop.pdf"),
                 width = plot_width, height = plot_height, limitsize = F)
      }else
        warning("No exclusive terms found for TW-based analysis with given filter")
      # -----

      # -----


    }
  }

}

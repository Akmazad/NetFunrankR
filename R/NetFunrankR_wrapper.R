#' Title This function works as a wrapper for running the whole functional analysis using network medicine based approach.
#'
#' @param disease_genes Input genes (a vector) that are used functional enrichment analysis. Default is "COPD-related genes".
#' @param PathwayRanker_Type A string that indicates the type of enrichment to be done. Options are: "ORA" (Hypergeometric-distribution based over-representation test), and "TW" (Topology-aware network medicine-based approach), "Both" (Default: both ORA and TW so that the results can be compared).
#' @param PE_ORA_pop A string that indicates the set of known annotaiton of pathways or functional term. Current version has the following options: "KEGG" (default), "Reactome" and "WikiPathway".
#' @param PE_ORA_cutoff_col A string to indicate the significance column for filter results based on a threshold. Options are: "p.adjust", "pvalue", or "NA" (default: no filter)
#' @param PE_ORA_cutoff_th A numeric value indicating the threshold (default: 0.05) for filtering result. Not effective if "PE_ORA_cutoff_col" parameter is set as "NA" (default)
#' @param nPermute A number (default: 1000) to indicate the number of permutation iteration to be conducted for hypothesis testing in Network-medicine based functional enrichment
#' @param neighbourhood_th A number (default: 1) to indicate the maximum value of shortest-path-length
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
    PathwayRanker_Type = "ORA", # "c('ORA','TW'")
    PE_ORA_pop = "KEGG", # c('KEGG', 'Reactome', 'WikiPathway')
    PE_ORA_cutoff_col = "p.adjust", # c("p.adjust", "pvalue")
    PE_ORA_cutoff_th = 0.05,
    nPermute = 10, # number of permutation iteration for hypothesis testing
    neighbourhood_th = 1, # maximum value of shortest-path-length
    string_PPI_score_th = 900, # controlling PPI confidence
    doPlot = TRUE,
    plot_width = 10,
    plot_height = 8,
    outdir = "./NetFunrankR_results/"
){
  if(is.null(disease_genes))
    stop("Disease genes are required for this analysis!!")

  set.seed(123)
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
    ## 2.1) do Pathway Enrichment test (ORA) [output-file]
    rankedPathways <- PathwayEnrichment_ORA(givenSet = COPD_genes,
                                            pathway.data = pathway.data,
                                            PE_ORA_pop = PE_ORA_pop)
    rankedPathways %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_ORA_result.csv"))
  }else{
    ## 2.2) Do Functional proximity [output-file]
    rankedPathways <- PathwayEnrichment_TopoAw(
      neighbourhood_th = neighbourhood_th,
      string_PPI_score_th = string_PPI_score_th,
      nPermute = nPermute,
      pathway.data = pathway.data,
      givenSet = COPD_genes)

    rankedPathways %>%
      fwrite(file = paste0(outdir, PE_ORA_pop, "_PE_TopoAw_result.csv"))
  }

}

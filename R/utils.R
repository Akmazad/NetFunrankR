PathwayEnrichment_ORA <- function(
    givenSet,
    PE_ORA_pop,
    pathway.data
    # cutoff_col = "pvalue",
    # cutoff_th = 0.05
){
  if(is.null(givenSet) ||
     length(givenSet) == 0 )
    stop("Pathway enrichment fails: it requires non-empty gene set")
  if(is.null(pathway.data) ||
     nrow(pathway.data) == 0)
    stop("Pathway enrichment fails: it requires Pathway annotation data")

  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(stringi)
  require(dplyr)


  # do the analysis -----
  pathwayNames = pathway.data$PathwayName %>% unique()
  nTerms <- pathwayNames %>% length()

  # background population is the whole human genome (protein-coding genes only)
  pop.genes <- get_all_human_protein_coding_genes()
  N = length(pop.genes)

  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(
    matrix(NA, nrow = nTerms, ncol = 7))

  K = length(givenSet)
  for (i in 1:nTerms) {
    # extract the genes in for a pathway
    pathway.genes.symbol = pathway.data %>%
      dplyr::filter(PathwayName == pathwayNames[i]) %>%
      dplyr::pull(GeneSymbol)

    # size of pathway gene set
    M = length(pathway.genes.symbol)

    # find the overlap
    x.overlap.genes.symbol = intersect(pathway.genes.symbol, givenSet)

    # size of the overlap
    x = length(x.overlap.genes.symbol)
    # if(x > 2){
      info.gene.overrep[i,1] = pathwayNames[i]
      info.gene.overrep[i,2] = paste0(pathway.genes.symbol, collapse = "/")
      info.gene.overrep[i,3] = paste0(x.overlap.genes.symbol, collapse = "/")

      x.overlap.genes.entrezID = clusterProfiler::bitr(x.overlap.genes.symbol,
                                                       fromType = "SYMBOL",
                                                       toType = "ENTREZID",
                                                       OrgDb = org.Hs.eg.db) %>%
        dplyr::select(ENTREZID) %>% unlist(use.names = F)


      info.gene.overrep[i,4] = paste0(x.overlap.genes.entrezID, collapse = "/")
      info.gene.overrep[i,5] = x
      info.gene.overrep[i,6] = x/M

      info.gene.overrep[i,7] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
    # }
  }
  # set the column names of the result
  colnames(info.gene.overrep) = c("Description",
                                  "Pathway.geneSymbol",
                                  "Overlapping.geneSymbol",
                                  "Overlapping.geneID",
                                  "Count", "GeneRatio",
                                  "pvalue")
  info.gene.overrep = info.gene.overrep %>%
    # dplyr::filter(pvalue != Inf) %>%
    dplyr::mutate("p.adjust" = p.adjust(pvalue,
                                        method = "fdr"
                                        # n = nrow(.)))
                                        ))
  res <- info.gene.overrep #%>%
    # dplyr::filter(!!dplyr::sym(cutoff_col) < cutoff_th) %>%
    # dplyr::arrange(!!dplyr::sym(cutoff_col))

  return(res)

}

PathwayEnrichment_TopoAw <- function(
    givenSet,
    pathway.data,
    neighbourhood_th,
    string_PPI_score_th,
    nPermute
    # cutoff_col = "p.adjust",
    # cutoff_th = 0.05
){
  if(is.null(pathway.data) ||
     nrow(pathway.data) == 0)
    stop("Pathway enrichment fails: it requires Pathway annotation data")

  # do the analysis -----
  pathwayNames = pathway.data$PathwayName %>% unique()
  nTerms <- pathwayNames %>% length()


  info.Pathway.network.proximity = data.frame(
    matrix(NA, nrow = nTerms, ncol = 4))

  for (i in 1:nTerms) {
  # for (i in 1:2) {
    pathwayname <- pathwayNames[i]

    # extract the genes in for a pathway
    pathway.genes.symbol = pathway.data %>%
      dplyr::filter(PathwayName == pathwayNames[i]) %>%
      dplyr::pull(GeneSymbol)

    # get network proximity score
    netP.scores <- quantify_network_proximity_score(
      neighbourhood_th = neighbourhood_th,
      string_PPI_score_th = string_PPI_score_th,
      module_a = givenSet,
      module_b = pathway.genes.symbol,
      nPermute = nPermute)

    # save the data
    info.Pathway.network.proximity[i,1] = pathwayname
    info.Pathway.network.proximity[i,2] = netP.scores[1] # observed distance
    info.Pathway.network.proximity[i,3] = netP.scores[2] # z-score
    info.Pathway.network.proximity[i,4] = netP.scores[3] # p-value

    cat(paste0("[",i,"] ", pathwayname, ": ",
               "; observed dis: ", netP.scores[1], " ",
               "; z-score: ", netP.scores[2], " ",
               "; pval: ", netP.scores[3], "\n"
    ))
  }
  # set the column names of the result
  colnames(info.Pathway.network.proximity) = c("Description",
                                  "observed_proximity",
                                  "z_score",
                                  "pvalue")
  info.Pathway.network.proximity = info.Pathway.network.proximity %>%
    # dplyr::filter(pvalue != Inf) %>%
    dplyr::mutate("p.adjust" = p.adjust(pvalue,
                                        method = "fdr"
                                        # n = nrow(.)))
                                        ))
  res <- info.Pathway.network.proximity #%>%
    # dplyr::filter(!!dplyr::sym(cutoff_col) < cutoff_th) %>%
    # dplyr::arrange(!!dplyr::sym(cutoff_col))

  return(res)
}

get_all_human_protein_coding_genes <- function(){
  require(EnsDb.Hsapiens.v75)

  edb <- EnsDb.Hsapiens.v75
  edb %>%
    AnnotationDbi::select(.,
                          keys=keys(edb, keytype="SYMBOL"),
                          columns=c("TXBIOTYPE"), keytype="SYMBOL") %>%
    dplyr::filter(TXBIOTYPE == "protein_coding") %>%
    dplyr::pull(SYMBOL) %>% unique() %>%
    return()
}

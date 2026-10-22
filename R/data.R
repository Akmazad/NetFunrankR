#' KEGG_2021_Human
#'
#' KEGG pathway annotation data from Enrichr Library
#'
#' @format A data frame with two variables:
#' \describe{
#' \item{\code{GeneSymbol}}{HUGO gene symbols belonging to a pathway}
#' \item{\code{PathwayName}}{Pathway name}
#' }
#' @source <https://maayanlab.cloud/Enrichr/#libraries>
"KEGG_2021_Human"

#' Reactome_Pathways_2024
#'
#' Reactome pathway annotation data from Enrichr Library
#'
#' @format A data frame with two variables:
#' \describe{
#' \item{\code{GeneSymbol}}{HUGO gene symbols belonging to a pathway}
#' \item{\code{PathwayName}}{Pathway name}
#' }
#' @source <https://maayanlab.cloud/Enrichr/#libraries>
"Reactome_Pathways_2024"

#' WikiPathways_2024_Human
#'
#' WikiPathway pathway annotation data from Enrichr Library
#'
#' @format A data frame with two variables:
#' \describe{
#' \item{\code{GeneSymbol}}{HUGO gene symbols belonging to a pathway}
#' \item{\code{PathwayName}}{Pathway name}
#' }
#' @source <https://maayanlab.cloud/Enrichr/#libraries>
"WikiPathways_2024_Human"

#' STRING_v12_PPI
#'
#' Human protein-protein associations from STRING v12.0 (taxon 9606), mapped
#' to HGNC symbols. Two scores are provided: the published combined score,
#' and a combined score recomputed (STRING channel-combination rule with
#' prior 0.041) from the neighborhood, fusion, co-occurrence, co-expression
#' and experimental channels only, i.e. without the curated-database and
#' text-mining channels that are partly derived from pathway databases.
#' Edges are retained if \code{nokb_score >= 400} or
#' \code{combined_score >= 700}.
#'
#' @format A data frame with four variables:
#' \describe{
#' \item{\code{from}, \code{to}}{HGNC gene symbols}
#' \item{\code{combined_score}}{STRING combined score (all channels)}
#' \item{\code{nokb_score}}{combined score without database/text-mining channels}
#' }
#' @source <https://string-db.org/> (protein.links.detailed.v12.0)
"STRING_v12_PPI"

#' COPD_proteinCodingGenes
#'
#' List of genes (HGNC symbols) found associated with COPD disease
#'
#' @format A character vector of HGNC gene symbols.
#' @source Morselli Gysi et al. (2021) <https://doi.org/10.1073/pnas.2025581118>
"COPD_proteinCodingGenes"

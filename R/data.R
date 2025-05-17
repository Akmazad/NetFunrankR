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
#' #' @source <https://maayanlab.cloud/Enrichr/#libraries>
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

#' protein.links
#'
#' Protein-Protein Interaction (human 9606) dataset from StringDB database
#'
#' @format A data frame with three variables:
#' \describe{
#' \item{\code{from}}{HUGO gene symbol 1}
#' \item{\code{to}}{HUGO gene symbol 2}
#' \item{\code{combined_score}}{Combined confidence score from multiple evidences}
#' }
#' @source <https://string-db.org/>
"protein.links"


#' COPD_proteinCodingGenes
#'
#' List of genes (HGNC symbols) found associated with COPD disease
#'
#' @format A vector:
#' \describe{
#' \item{\code{value}}{Hugo Gene symbol}
#' }
#' @source <https://doi.org/10.1073/pnas.2301342120>
"COPD_proteinCodingGenes"

#' KEGG pathways (Enrichr KEGG_2021_Human)
#'
#' @format A \code{data.frame} with columns \code{GeneSymbol} (HGNC symbol)
#'   and \code{PathwayName}; 320 pathways.
#' @source \url{https://maayanlab.cloud/Enrichr/#libraries}; see
#'   \code{inst/script/make-STRING_v12_PPI.R}.
#' @seealso \code{\link{loadPathways}}
#' @examples
#' data(KEGG_2021_Human, package = "NetFunrankR")
#' length(unique(KEGG_2021_Human$PathwayName))
"KEGG_2021_Human"

#' Reactome pathways (Enrichr Reactome_Pathways_2024)
#'
#' @format A \code{data.frame} with columns \code{GeneSymbol} and
#'   \code{PathwayName}.
#' @source \url{https://maayanlab.cloud/Enrichr/#libraries}
#' @seealso \code{\link{loadPathways}}
#' @examples
#' data(Reactome_Pathways_2024, package = "NetFunrankR")
#' length(unique(Reactome_Pathways_2024$PathwayName))
"Reactome_Pathways_2024"

#' WikiPathways (Enrichr WikiPathways_2024_Human)
#'
#' @format A \code{data.frame} with columns \code{GeneSymbol} and
#'   \code{PathwayName}.
#' @source \url{https://maayanlab.cloud/Enrichr/#libraries}
#' @seealso \code{\link{loadPathways}}
#' @examples
#' data(WikiPathways_2024_Human, package = "NetFunrankR")
#' length(unique(WikiPathways_2024_Human$PathwayName))
"WikiPathways_2024_Human"

#' STRING v12 human protein associations
#'
#' Human protein associations from STRING v12.0 mapped to HGNC symbols, with
#' two scores: the published combined score, and a combined score recomputed
#' with STRING's channel-combination rule (prior 0.041) from the
#' neighbourhood, fusion, co-occurrence, co-expression and experimental
#' channels only, i.e. without the curated-database and text-mining channels
#' that are partly derived from pathway databases. Edges are retained if
#' \code{nokb_score >= 400} or \code{combined_score >= 700}.
#'
#' @format A \code{data.frame} with columns \code{from}, \code{to} (HGNC
#'   symbols), \code{combined_score} and \code{nokb_score} (0--1000).
#' @source STRING v12.0 (\url{https://string-db.org}), CC BY 4.0; created
#'   with \code{inst/script/make-STRING_v12_PPI.R}.
#' @references Szklarczyk D, et al. (2023) The STRING database in 2023.
#'   \emph{Nucleic Acids Research} 51:D638--D646.
#' @seealso \code{\link{stringNetwork}}
#' @examples
#' data(STRING_v12_PPI, package = "NetFunrankR")
#' summary(STRING_v12_PPI$nokb_score)
"STRING_v12_PPI"

#' COPD-associated protein-coding genes
#'
#' @format A character vector of 307 HGNC gene symbols.
#' @source Morselli Gysi D, et al. (2021) Network medicine framework for
#'   identifying drug-repurposing opportunities for COVID-19.
#'   \emph{PNAS} 118:e2025581118.
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' head(COPD_proteinCodingGenes)
"COPD_proteinCodingGenes"

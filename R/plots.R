## Colours of the evidence classes (validated categorical palette; "none"
## is neutral grey). Exported via evidenceColors() for reuse.
.evidenceCols <- c("overlap + network" = "#1baf7a",
    "overlap only" = "#2a78d6", "network only" = "#eb6834",
    "joint" = "#6250d6", "none" = "#b9b8b2")

#' Colours used for the evidence classes
#'
#' @returns named character vector of hex colours.
#' @export
#' @examples
#' evidenceColors()
evidenceColors <- function() .evidenceCols

#' Evidence map: overlap versus beyond-overlap network evidence
#'
#' Places each gene set by the significance of its degree-aware overlap test
#' (x) and of its beyond-overlap test (y), so that sets supported by direct
#' membership, by network proximity of non-member query genes, or by both
#' fall in different regions.
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param alpha significance threshold (evidence classes).
#' @param nLabel number of top gene sets (by \code{p_combined}) to label.
#' @returns a \code{ggplot} object.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' plotEvidenceMap(res)
plotEvidenceMap <- function(x, alpha = 0.05, nLabel = 12L) {
    df <- as.data.frame(x, alpha = alpha)
    yCol <- if (x@params$beyond == "distance") "p_beyond" else "p_linker"
    df$xv <- -log10(df$p_overlap)
    df$yv <- -log10(df[[yCol]])
    df$evidence <- factor(df$evidence, .evidenceLevels)
    lab <- utils::head(df[order(df$p_combined), ], nLabel)
    ggplot2::ggplot(df, ggplot2::aes(.data$xv, .data$yv)) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$evidence,
            size = .data$set_in_net), alpha = 0.8) +
        ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$set),
            size = 2.8, vjust = -0.9, check_overlap = TRUE) +
        ggplot2::scale_colour_manual(values = .evidenceCols, drop = TRUE,
            name = NULL) +
        ggplot2::scale_size_area(max_size = 6, name = "Set size") +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.03, 0.25))) +
        ggplot2::labs(
            x = expression(-log[10] ~ P[overlap] ~ "(degree-aware)"),
            y = if (yCol == "p_linker")
                expression(-log[10] ~ P[linker] ~ "(beyond overlap)") else
                expression(-log[10] ~ P[distance] ~ "(beyond overlap)")) +
        ggplot2::theme_minimal(base_size = 11)
}

#' Dot plot of the top gene sets
#'
#' @inheritParams plotEvidenceMap
#' @param n number of gene sets to show.
#' @returns a \code{ggplot} object.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork())
#' plotTopSets(res, n = 15)
plotTopSets <- function(x, n = 20L, alpha = 0.05) {
    df <- utils::head(as.data.frame(x, alpha = alpha), n)
    df$set <- factor(df$set, rev(df$set))
    df$evidence <- factor(df$evidence, .evidenceLevels)
    ggplot2::ggplot(df, ggplot2::aes(-log10(.data$p_combined), .data$set)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend =
            -log10(.data$p_combined), yend = .data$set), colour = "grey80") +
        ggplot2::geom_point(ggplot2::aes(colour = .data$evidence,
            size = .data$overlap + .data$n_linker)) +
        ggplot2::scale_colour_manual(values = .evidenceCols, drop = TRUE,
            name = NULL) +
        ggplot2::scale_size_area(max_size = 6,
            name = "Overlap + linker genes") +
        ggplot2::labs(x = expression(-log[10] ~ P[combined]), y = NULL) +
        ggplot2::theme_minimal(base_size = 11)
}

#' Exact null distribution of the proximity statistic
#'
#' @param nd output of \code{\link{nullDistribution}}.
#' @param which \code{"S"} (total proximity) or \code{"Sb"} (beyond-overlap
#'   distance).
#' @param title plot title.
#' @returns a \code{ggplot} object.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' nd <- nullDistribution(COPD_proteinCodingGenes, loadPathways("KEGG"),
#'     stringNetwork(), "Nicotine addiction")
#' plotNullDistribution(nd)
plotNullDistribution <- function(nd, which = c("S", "Sb"), title = NULL) {
    which <- match.arg(which)
    p <- nd[[which]]
    denom <- if (which == "S") nd$n else nd$n - nd$overlap
    df <- data.frame(d = (seq_along(p) - 1) / denom, prob = p)
    df <- df[df$prob > 1e-8 * max(df$prob), ]
    ggplot2::ggplot(df, ggplot2::aes(.data$d, .data$prob)) +
        ggplot2::geom_col(fill = "grey65", width = 0.8 / denom) +
        ggplot2::geom_vline(xintercept = nd$S_obs / denom,
            colour = "#eb6834", linewidth = 1) +
        ggplot2::labs(x = "Mean shortest-path distance to the gene set",
            y = "Null probability", title = title) +
        ggplot2::theme_minimal(base_size = 11)
}

#' Network view of query genes linked to a gene set
#'
#' Draws the query genes that are members of, or adjacent to, a gene set
#' together with the set members they interact with.
#'
#' @param x an \code{\linkS4class{NFResult}}.
#' @param set name of a tested gene set.
#' @param network the \code{\linkS4class{NFNetwork}} used for \code{x}.
#' @param maxMembers maximum number of set members to draw (most connected
#'   to the query first).
#' @returns a \code{ggplot} object.
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"), net)
#' plotLinkerNetwork(res, "Nicotine addiction", net)
plotLinkerNetwork <- function(x, set, network, maxMembers = 40L) {
    sub <- linkerSubnetwork(x, set, network, maxMembers)
    lay <- .stableLayout(sub$graph)
    nd <- data.frame(name = igraph::V(sub$graph)$name, x = lay[, 1],
        y = lay[, 2])
    nd$role <- sub$role[nd$name]
    el <- igraph::as_edgelist(sub$graph)
    ed <- data.frame(x = nd$x[match(el[, 1], nd$name)],
        y = nd$y[match(el[, 1], nd$name)],
        xend = nd$x[match(el[, 2], nd$name)],
        yend = nd$y[match(el[, 2], nd$name)])
    cols <- c("Query gene in set" = "#1baf7a",
        "Query linker gene" = "#eb6834", "Set member" = "#9fb6c3")
    ggplot2::ggplot() +
        ggplot2::geom_segment(data = ed, ggplot2::aes(.data$x, .data$y,
            xend = .data$xend, yend = .data$yend), colour = "grey80",
            linewidth = 0.3) +
        ggplot2::geom_point(data = nd, ggplot2::aes(.data$x, .data$y,
            fill = .data$role), shape = 21, size = 4, colour = "white") +
        ggplot2::geom_text(data = nd[nd$role != "Set member", ],
            ggplot2::aes(.data$x, .data$y, label = .data$name),
            size = 2.6, vjust = -1.1) +
        ggplot2::scale_fill_manual(values = cols, name = NULL) +
        ggplot2::labs(title = set) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(legend.position = "bottom")
}

## Force-directed (graphopt) layout started from circle coordinates:
## deterministic and does not use (or change) the user's random-number state.
.stableLayout <- function(g) {
    igraph::layout_with_graphopt(g, start = igraph::layout_in_circle(g))
}

#' Subnetwork of query genes linked to a gene set
#'
#' @inheritParams plotLinkerNetwork
#' @returns list with \code{graph} (an \code{igraph} induced subgraph) and
#'   \code{role} (named character vector: "Query gene in set",
#'   "Query linker gene" or "Set member").
#' @export
#' @examples
#' data(COPD_proteinCodingGenes, package = "NetFunrankR")
#' net <- stringNetwork()
#' res <- netEnrich(COPD_proteinCodingGenes, loadPathways("KEGG"), net)
#' sub <- linkerSubnetwork(res, "Nicotine addiction", net)
#' table(sub$role)
linkerSubnetwork <- function(x, set, network, maxMembers = 40L) {
    Y <- x@geneSets[[set]]
    if (is.null(Y))
        stop("gene set '", set, "' was not tested in this result")
    g <- nfGraph(network)
    vn <- igraph::V(g)$name
    Y <- intersect(Y, vn)
    q <- x@query
    nb <- igraph::adjacent_vertices(g, setdiff(q, Y))
    hits <- lapply(nb, function(v) intersect(vn[v], Y))
    linkers <- names(hits)[lengths(hits) > 0]
    score <- sort(table(unlist(hits[linkers])), decreasing = TRUE)
    inSet <- intersect(q, Y)
    ## every linker keeps its best-connected partner; remaining slots go to
    ## the members most connected to the query
    partner <- vapply(hits[linkers], function(h)
        h[which.max(score[h])], character(1))
    members <- unique(c(setdiff(partner, inSet),
        setdiff(names(score), inSet)))
    members <- utils::head(members, max(maxMembers, length(unique(partner))))
    nodes <- unique(c(inSet, linkers, members))
    role <- stats::setNames(rep("Set member", length(nodes)), nodes)
    role[linkers] <- "Query linker gene"
    role[inSet] <- "Query gene in set"
    list(graph = igraph::induced_subgraph(g, nodes), role = role)
}

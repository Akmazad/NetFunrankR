.nf_colors <- c("overlap + network" = "#7A4E9E",
                "overlap only" = "#616530",
                "network only" = "#CC8214",
                "joint" = "#3B7EA1",
                "none" = "grey75")

#' Venn diagram of gene sets significant by two methods
#' @param vec1,vec2 character vectors of significant set names.
#' @param names labels of the two methods.
#' @returns a ggplot object.
#' @export
plot_venn <- function(vec1, vec2,
                      names = c("Over-representation", "Network Proximity")) {
  x <- list(vec1, vec2)
  names(x) <- names
  ggvenn::ggvenn(x, fill_color = c("#616530", "#CC8214"),
                 stroke_size = 0.1, set_name_size = 5, text_size = 5,
                 stroke_alpha = 0.8)
}

#' Evidence map: overlap versus beyond-overlap network evidence
#'
#' Each gene set is placed by the significance of its degree-aware overlap
#' test (x) and of its conditional beyond-overlap proximity test (y), so that
#' sets supported by direct membership, by network proximity of non-member
#' query genes, or by both, fall in different quadrants.
#'
#' @param res an \code{nf_result}.
#' @param alpha FDR threshold used for the quadrant lines and colours.
#' @param n_label number of top sets (by \code{p_combined}) to label.
#' @returns a ggplot object.
#' @export
plot_evidence_map <- function(res, alpha = 0.05, n_label = 12) {
  df <- as.data.frame(res)
  df$class <- nf_evidence_class(res, alpha)
  df$x <- -log10(df$p_overlap)
  df$y <- -log10(df$p_linker)
  lab <- df[order(df$p_combined), ][seq_len(min(n_label, nrow(df))), ]
  thr <- function(p, fdr) {
    s <- p[fdr < alpha]
    if (length(s)) -log10(max(s)) else NA_real_
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_vline(xintercept = thr(df$p_overlap, df$fdr_overlap),
                        linetype = "dashed", colour = "grey50", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = thr(df$p_linker, df$fdr_linker),
                        linetype = "dashed", colour = "grey50", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$class, size = .data$set_in_net),
                        alpha = 0.8) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$set),
                       size = 2.8, vjust = -0.9, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = .nf_colors, drop = TRUE, name = NULL) +
    ggplot2::scale_size_area(max_size = 6, name = "Set size") +
    ggplot2::labs(x = expression(-log[10]~P[overlap]~"(degree-aware)"),
                  y = expression(-log[10]~P[linker]~"(beyond overlap)")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.25))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right")
}

#' Exact null distribution of the proximity statistic for one gene set
#'
#' @param nd output of \code{\link{nf_null_distribution}}.
#' @param which "S" (total proximity) or "Sb" (beyond-overlap proximity).
#' @param title plot title.
#' @returns a ggplot object.
#' @export
plot_null <- function(nd, which = c("S", "Sb"), title = NULL) {
  which <- match.arg(which)
  p <- nd[[which]]
  denom <- if (which == "S") nd$n else nd$n - nd$overlap
  df <- data.frame(d = (seq_along(p) - 1) / denom, prob = p)
  df <- df[df$prob > 1e-8 * max(df$prob), ]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$d, y = .data$prob)) +
    ggplot2::geom_col(fill = "grey60", width = 0.8 / denom) +
    ggplot2::geom_vline(xintercept = nd$S_obs / denom, colour = "#CC8214",
                        linewidth = 1) +
    ggplot2::labs(x = "Mean shortest-path distance to gene set",
                  y = "Null probability", title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Network view of query genes linked to a gene set
#'
#' Draws the query genes that are members of, or adjacent to, the gene set,
#' together with the set members they touch.
#'
#' @param query query genes.
#' @param set_genes genes of the gene set.
#' @param net an \code{nf_network}.
#' @param max_members maximum number of set members to draw (highest number
#'   of links to the query first).
#' @param seed layout seed.
#' @returns a ggplot object.
#' @export
plot_linker_network <- function(query, set_genes, net, max_members = 40, seed = 1) {
  g <- net$graph
  q <- intersect(query, net$genes)
  Y <- intersect(set_genes, net$genes)
  qn <- setdiff(q, Y)
  nb <- igraph::adjacent_vertices(g, qn)
  linkers <- qn[vapply(nb, function(v) any(igraph::V(g)$name[v] %in% Y), logical(1))]
  touched <- unique(unlist(lapply(nb[linkers], function(v) intersect(igraph::V(g)$name[v], Y))))
  members_q <- intersect(q, Y)
  score <- table(factor(unlist(lapply(nb[linkers], function(v)
    intersect(igraph::V(g)$name[v], Y))), levels = touched))
  keepY <- union(members_q, names(sort(score, decreasing = TRUE))[seq_len(min(max_members, length(score)))])
  nodes <- unique(c(linkers, keepY))
  sg <- igraph::induced_subgraph(g, nodes)
  set.seed(seed)
  lay <- igraph::layout_with_fr(sg)
  vn <- igraph::V(sg)$name
  role <- ifelse(vn %in% members_q, "Query gene in set",
                 ifelse(vn %in% linkers, "Query linker gene", "Set member"))
  nd <- data.frame(name = vn, x = lay[, 1], y = lay[, 2], role = role)
  el <- igraph::as_edgelist(sg)
  ed <- data.frame(x = nd$x[match(el[, 1], vn)], y = nd$y[match(el[, 1], vn)],
                   xend = nd$x[match(el[, 2], vn)], yend = nd$y[match(el[, 2], vn)])
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed, ggplot2::aes(x = .data$x, y = .data$y,
                                                  xend = .data$xend, yend = .data$yend),
                          colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_point(data = nd, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$role),
                        shape = 21, size = 4, colour = "white") +
    ggplot2::geom_text(data = nd[nd$role != "Set member", ],
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
                       size = 2.6, vjust = -1.1) +
    ggplot2::scale_fill_manual(values = c("Query gene in set" = "#7A4E9E",
                                          "Query linker gene" = "#CC8214",
                                          "Set member" = "#9DB4C0"), name = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

#' Mirrored bar plot of p-values for gene sets found by two methods
#'
#' @param df data.frame with columns \code{set}, \code{p1}, \code{p2}.
#' @param names labels of the two methods.
#' @returns a ggplot object.
#' @export
plot_double_bar <- function(df, names = c("Over-representation", "Network Proximity")) {
  df$sort_key <- pmin(df$p1, df$p2)
  long <- rbind(data.frame(set = df$set, group = names[1], logp = log10(df$p1), key = df$sort_key),
                data.frame(set = df$set, group = names[2], logp = -log10(df$p2), key = df$sort_key))
  long$set <- stats::reorder(long$set, -long$key)
  m <- max(abs(long$logp), na.rm = TRUE)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$logp, y = .data$set, fill = .data$group)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = c(log10(0.05), -log10(0.05)),
                        linetype = "dashed", colour = "red") +
    ggplot2::scale_x_continuous(limits = c(-m, m), labels = function(x) abs(x)) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#616530", "#CC8214"), names)) +
    ggplot2::labs(x = expression(-log[10]*"(P-value)"), y = NULL, fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", panel.grid.major.y = ggplot2::element_blank())
}

#' Lollipop plot of gene-set p-values
#'
#' @param df data.frame with columns \code{set} and \code{pvalue}.
#' @param color_palette two colours (high, low).
#' @returns a ggplot object.
#' @export
plot_lollipop <- function(df, color_palette = c("#CC8214", "#f6e8b1")) {
  df$logp <- -log10(df$pvalue)
  df$set <- stats::reorder(df$set, df$logp)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$logp, y = .data$set)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$logp, yend = .data$set),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$logp), size = 4) +
    ggplot2::scale_colour_gradient(high = color_palette[1], low = color_palette[2],
                                   name = expression(-log[10](P))) +
    ggplot2::labs(x = expression(-log[10](P-value)), y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

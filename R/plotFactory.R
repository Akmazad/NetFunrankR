plot_venn <- function(vec1, vec2){
  # require(ggvenn)

  x <- list(
    A = vec1,
    B = vec2
  )
  names(x) <- c("Over-representation","Network Proximity")
  ggvenn(
    x,
    fill_color = c("#616530FF", "#CC8214FF"),
    stroke_size = 0.1, set_name_size = 7.5, padding = 0.05, text_size = 7.5, stroke_alpha = 0.8
  )
}

plot_double_bar <- function(df) {
  # # Required packages
  # if (!requireNamespace("ggplot2", quietly = TRUE) ||
  #     !requireNamespace("tidyr", quietly = TRUE) ||
  #     !requireNamespace("forcats", quietly = TRUE) ||
  #     !requireNamespace("dplyr", quietly = TRUE)) {
  #   stop("Please install ggplot2, tidyr, forcats, and dplyr packages.")
  # }
  #
  # library(dplyr)
  # library(tidyr)
  # library(forcats)
  # library(ggplot2)

  # Threshold for significance
  log_threshold <- -log10(0.05)

  # Add sort_key as the minimum of pvalue.x and pvalue.y
  df <- df %>%
    mutate(sort_key = pmin(pvalue.x, pvalue.y, na.rm = TRUE))

  # Reshape and compute mirrored -log10(pvalue)
  df_long <- df %>%
    pivot_longer(
      cols = c(pvalue.x, pvalue.y),
      names_to = "group",
      values_to = "pvalue"
    ) %>%
    mutate(
      logp = -log10(pvalue),
      logp = ifelse(group == "pvalue.x", -logp, logp),
      group = recode(group,
                     "pvalue.x" = "Over-representation",
                     "pvalue.y" = "Network Proximity"),
      Description = fct_reorder(Description, df$sort_key[match(Description, df$Description)])
    )

  # compute symmetric x-limits
  max_abs <- max(abs(df_long$logp), na.rm = TRUE)

  # Plot
  ggplot(df_long, aes(x = logp, y = Description, fill = group)) +
    geom_col(width = 0.7) +
    # vertical thresholds
    geom_vline(xintercept = c(-log_threshold, log_threshold),
               linetype = "dashed", color = "red") +
    # symmetric scale, labels positive
    scale_x_continuous(
      limits = c(-max_abs, max_abs),
      labels = function(x) abs(x),
      expand = expansion(add = 0)
    ) +
    scale_fill_manual(
      values = c("Over-representation" = "#616530FF",
                 "Network Proximity" = "#CC8214FF")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = expression(-log[10]*"(P-value)"),
      y = NULL,
      fill = NULL,
      title = ""
    ) +
    theme_minimal() +
    theme(
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.title       = element_blank(),
      legend.justification = "center",
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

plot_lollipop <- function(df, color_palette) {
  # # check for required packages
  # if (!requireNamespace("ggplot2", quietly = TRUE) ||
  #     !requireNamespace("dplyr", quietly = TRUE) ||
  #     !requireNamespace("forcats", quietly = TRUE)) {
  #   stop("Please install ggplot2, dplyr and forcats.")
  # }
  #
  # library(dplyr)
  # library(forcats)
  # library(ggplot2)

  # prepare the data
  df2 <- df %>%
    mutate(
      logp = -log10(pvalue),
      Description = fct_reorder(Description, logp)
    )

  # find min/max for the color scale
  logp_range <- range(df2$logp, na.rm = TRUE)

  # plot
  ggplot(df2, aes(x = logp, y = Description)) +
    # stems
    geom_segment(aes(x = 0, xend = logp, y = Description, yend = Description),
                 color = "grey70") +
    # lollipop heads, colored by logp
    geom_point(aes(color = logp), size = 4) +
    # continuous gradient from low to high
    scale_color_gradient(
      high = color_palette[1],
      low  = color_palette[2],
      limits = logp_range,
      name = expression(-log[10](P-value))
    ) +
    labs(
      x = expression(-log[10](P-value)),
      y = NULL,
      title = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major.y = element_blank()
    )
}


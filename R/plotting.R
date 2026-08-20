boxplot_with_pairwise_pvalues <- function(value, group, comparisons, colors,
                                           test = "t.test", title = NULL,
                                           ylab = NULL, xlab = NULL,
                                           rotate_x_labels = TRUE,
                                           add_points = TRUE) {
  library(ggplot2)
  library(ggpubr)

  df <- data.frame(group = group, value = as.numeric(value))
  df <- df[order(df$group), ]

  p <- ggboxplot(df, x = "group", y = "value", fill = "group", palette = colors,
                  add = if (add_points) "jitter" else "none") +
    stat_compare_means(method = test, comparisons = comparisons, label = "p.format") +
    labs(title = title, x = xlab, y = ylab) +
    theme(legend.position = "none")

  if (rotate_x_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  }
  p
}

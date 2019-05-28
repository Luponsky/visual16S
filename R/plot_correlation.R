#' plot_correlation
#'
#' plot_correlation can plot correlation using a correlation table.
#'
#' @param cor_tab A data frame that contains at least two columns to calculate
#' correlation.
#'
#' @param x Colname of cor_tab. x and y must have the same length.
#' @param y Colname of cor_tab. x and y must have the same length.
#'
#' @param method A character string indicating which correlation coefficient is
#' to be used. One of c("pearson", "kendall", or "spearman"), default is
#' "pearson".
#'
#' @export

plot_correlation <- function (cor_tab, x, y, method = "pearson") {
  # Notice: Colnames of the input table can only be letters or numbers.
  if (any(str_detect(c(x, y), '\\W'))) {
    stop(paste0("Colnames of the input columns can only contain letters or",
                " numbers, or it can't be recognized when plotting."))
  }
  if (method == "pearson") {
    unit <- "Pearson's r"
  } else if (method == "spearman") {
    unit <- "Spearman's Rho"
  } else if (method == "kendall") {
    unit <- "Kendall's Tau"
  } else {
    stop("Input method is not supported.")
  }
  # Calculate correlation
  cor_res <- cor.test(cor_tab[[x]], cor_tab[[y]], method = method)
  if (cor_res$p.value < 2.2e-16) {
    ggplot(data = cor_tab, aes_string(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = lm) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.1,
               label = paste0(unit, " = ", round(cor_res$estimate, 2))) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.05,
               label = paste0('p-value < 2.2e-16')) +
      labs(x = as.character(x), y = as.character(y)) +
      theme_bw() +
      theme(panel.grid = element_blank())
  } else {
    ggplot(data = cor_tab, aes_string(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = lm) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.1,
               label = paste0(unit, " = ", round(cor_res$estimate, 2))) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.05,
               label = paste0('p-value = ', round(cor_res$p.value, 3))) +
      labs(x = as.character(x), y = as.character(y)) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
}
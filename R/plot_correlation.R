#' plot_correlation
#'
#' plot_correlation can plot correlation using a correlation table.
#'
#' @param cor_tab A data frame that contains at least two columns to calculate correlation.
#'
#' @param x Colname of cor_tab. If x's length is more than 1, plot_correlation will calculate
#' correlation respectively and facet the plot by x.
#'
#' @param y Colname of cor_tab.
#'
#' @param method A character string indicating which correlation coefficient is to be used. One of
#' c("pearson", "kendall", or "spearman"), default is "pearson".
#'
#' @export

plot_correlation <- function (cor_tab, x, y, method = "pearson") {
  # Notice: Colnames of the input table can only be letters or numbers.
  #if (any(str_detect(c(x, y), '\\W'))) {
  #  stop(paste0("Colnames of the input columns can only contain letters or numbers, or it can't be ",
  #              "recognized when plotting."))
  #}
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
  if (length(x) == 1) {
    cor_res <- cor.test(cor_tab[[x]], cor_tab[[y]], method = method)
    if (cor_res$p.value < 2.2e-16) {
      p_value <- "p-value < 2.2e-16"
    } else {
      if (cor_res$p.value < 0.001) {
        p_value <- paste0("p-value = ", formatC(cor_res$p.value, format = "e", digits = 1))
      } else {
        p_value <- paste0('p-value = ', round(cor_res$p.value, 3))
      }
    }
    g <- ggplot(data = cor_tab, aes_string(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = lm) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.1,
               label = paste0(unit, " = ", round(cor_res$estimate, 2))) +
      annotate(geom = 'text',
               x = max(cor_tab[[x]]) / 2,
               y = max(cor_tab[[y]]) * 1.05,
               label = p_value) +
      labs(x = as.character(x), y = as.character(y)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12))
    ggExtra::ggMarginal(g, type = "histogram", fill = "transparent")
  } else {
    cor_res <- data.frame(row.names = c("facet", "correlation", "pvalue", "x")) %>%
      t() %>% as.data.frame()
    for (i in 1:length(x)) {
      z <- cor.test(cor_tab[[x[i]]], cor_tab[[y]], method = method)
      cor_res[i,1] <- x[i]
      cor_res[i,2] <- paste0(unit, " = ", round(as.numeric(z$estimate), 2))
      if (z$p.value < 2.2e-16) {
        cor_res[i,3] <- "p-value < 2.2e-16"
      } else if (z$p.value < 0.001) {
        cor_res[i,3] <- paste0("p-value = ", formatC(z$p.value, format = "e", digits = 1))
      } else {
        cor_res[i,3] <- paste0('p-value = ', round(z$p.value, 3))
      }
      cor_res[i,4] <- max(cor_tab[[x[i]]]) / 2
    }
    cor_tab <- gather(cor_tab, x, key = "facet", value = "values")
    values <- "values"
    ggplot(cor_tab, aes_string(x = values, y = y)) +
      facet_wrap(vars(facet), scales = "free") +
      geom_point() +
      geom_smooth(method = lm) +
      xlab("") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12)) +
      geom_text(data = cor_res, mapping = aes(x = x, y = Inf, label = correlation),
                vjust = 2) +
      geom_text(data = cor_res, mapping = aes(x = x, y = Inf, label = pvalue),
                vjust = 4)
  }
}

#' plot_alpha_diversity
#'
#' This is a function for plotting alpha diversity.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#' metadata and phylogenetic tree.
#' @param feature The column name of the feature you want to select from
#' metadata.
#' @param feature2 The column name of another feature you want to select from
#' metadata, e.g. "Gender", which will make the plots draw in different shapes.
#' Default is NA.
#' @param measures The measures to calculate alpha diversity. Default is NA. If
#' NA, all available alpha diversity measures will be calculated and generate a
#' table. If not NA, measures should be one of c("Observed", "Chao1", "ACE",
#' "Shannon", "Simpson", "InvSimpson", "Fisher").
#' @param p_test The p-value to test alpha diversity. p_test should be either
#' "wilcox" or "kruskal". PS: "wilcox" can only work with two groups.
#' @export
#' @examples
#' plot_alpha_diversity(demo_phyloseq_object, feature = "diagnosis",
#'                      feature2 = NA, measures = "Chao1", p_test = "kruskal")

plot_alpha_diversity <- function (phyloseq, feature, feature2 = NA,
                                  measures = NA, p_test = "kruskal") {
  set.seed(99)
  if (!is.na(measures)) {
    if (!measures %in% c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                         "InvSimpson", "Fisher")) {
      stop(paste0('Argument "measures" should be one of c("Observed", "Chao1"',
                  ', "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher").'))
    } else {
      ## Step 1: Use plot_richness function to calculate alpha diversity
      alpha_diversity <- plot_richness(phyloseq, x = feature,
                                       measures = measures)
      ## Step 2: Calculate p-value
      if (p_test == "wilcox") {
        # Prepare feature table for calculating Mann-Whitney U test
        feature_tab_4_MWtest <- extract_metadata_phyloseq(phyloseq, feature)
        # Extract feature levels
        feature_0 <- feature_tab_4_MWtest[[feature]] %>% unique() %>% .[1]
        feature_1 <- feature_tab_4_MWtest[[feature]] %>% unique() %>% .[2]
        replace_feature <- c(as.character(feature_0), as.character(feature_1))
        # Revalue feature levels as 0 and 1
        feature_tab_4_MWtest[[feature]] <-
          plyr::mapvalues(feature_tab_4_MWtest[[feature]], to = c(0, 1),
                          from = replace_feature)
        p_value <- wilcox.test(alpha_diversity$data$value ~
                                 feature_tab_4_MWtest[[feature]])$p.value
      } else if (p_test == "kruskal") {
        # Kruskal test
        p_value <- kruskal.test(alpha_diversity$data$value,
                                factor(alpha_diversity$data[,feature]))$p.value
      } else {
        stop("The input p_test is not supported.")
      }
      ## Step 3: Plot alpha diversity
      y <- "value"
      if (is.na(feature2)) {
        p <- ggplot(data = alpha_diversity$data,
                    # Use aes_string() to pass variables to ggplot
                    aes_string(x = feature, y = y, color = feature)) +
          geom_boxplot() +
          geom_jitter(size = 3) +
          ylab(paste0(measures, " Diversity")) +
          annotate("text",
                   x = ((alpha_diversity$data[[feature]] %>% unique() %>%
                           length() + 1)/2),
                   y = (max(alpha_diversity$data$value) * 1.1),
                   label = paste0(p_test, " p-value = ", round(p_value, 3)),
                   size = 3) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                legend.text = element_text(size = 12),
                strip.text.x = element_text(size = 14))
        p + ggsci::scale_color_jco() + ggsci::scale_fill_jco()
      } else {
        p <- ggplot(data = alpha_diversity$data,
                    # Use aes_string() to pass variables to ggplot
                    aes_string(x = feature, y = y, color = feature)) +
          geom_boxplot() +
          geom_jitter(size = 3, aes_string(shape = feature2)) +
          ylab(paste0(measures, " Diversity")) +
          annotate("text",
                   x = ((alpha_diversity$data[[feature]] %>% unique() %>%
                           length() + 1)/2),
                   y = (max(alpha_diversity$data$value) * 1.1),
                   label = paste0(p_test, " p-value = ", round(p_value, 3)),
                   size = 3) +
          scale_shape_manual(values = c(0:6)) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                legend.text = element_text(size = 12),
                strip.text.x = element_text(size = 14))
        p + ggsci::scale_color_jco() + ggsci::scale_fill_jco()
      }
    }
  } else {
    alpha_diversity <- plot_richness(phyloseq, x = feature) %>%
      .[["data"]] %>%
      select(-se) %>%
      spread(variable, value)
    return(alpha_diversity)
  }
}

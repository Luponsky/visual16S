#' plot_stacked_bar
#'
#' This is a function for stacked bar plot.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#'                 metadata and phylogenetic tree.
#' @param feature The feature that shows in x-axis text with different colors.
#' @param level Which taxonomy level to plot.
#' @param x_size The front size of x-axis text. Default is 8.
#' @param legend_position Legend position. Default is "top". legend_position
#'                        should be one of "none", "left", "right", "bottom",
#'                        "top".
#' @param legend_size Lengend size. Default is 10.
#' @export
#' @examples
#' plot_stacked_bar(demo_phyloseq_object, feature = "diagnosis",
#'                  level = "Order")

plot_stacked_bar <- function(phyloseq, feature, level, x_size = 8,
                             legend_position = "top", legend_size = 10) {
  ## Step 1: First construct otu then convert to percentage
  otu_percent <- construct_otu_table(phyloseq, level) %>%
    convert_to_percentage(row_sum = TRUE) %>%
    as.data.frame()
  ## Step 2: Construct table for stacked bar plot
  plot_tab <- otu_percent %>%
    # First re-order taxonomy by total counts
    .[,order(colSums(.), decreasing = TRUE)] %>%
    # Then re-order sample by the most abundance taxonomy
    .[order(.[,1], decreasing = TRUE),] %>%
    # Turn SampleID to a column
    rownames_to_column(var = "SampleID")
  # Prepare levels for SampleID
  levels_SampleID <- plot_tab$SampleID
  # Prepare levels for taxonomy levels
  levels_level <- colnames(plot_tab)[2:ncol(plot_tab)]
  ## Step 3: Join metadata
  sample_feature <- extract_metadata_phyloseq(phyloseq, feature)
  plot_tab <- left_join(plot_tab, sample_feature)
  # Prepare levels for feature
  levels_feature <- plot_tab[,ncol(plot_tab)]
  # Add levels to SampleID
  plot_tab$SampleID <- factor(plot_tab$SampleID, levels = levels_SampleID)
  ## Step 4: Turn plot_table to a long table for plotting
  plot_tab <- gather(plot_tab,
                     colnames(plot_tab)[2:(ncol(plot_tab)-1)],
                     key = level,
                     value = "abundance")
  # Add levels to taxonomy levels
  plot_tab$level <- factor(plot_tab$level, levels = levels_level)
  ## Step 5: Bar plot
    ggplot(plot_tab, aes(x = SampleID, y = abundance)) +
      scale_fill_manual(values = distinctive_colors) +
      geom_bar(mapping = aes(fill = level), position = "fill",
               stat = "identity") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            axis.text.x = element_text(size = x_size, angle = 90,
                                       color = levels_feature),
            legend.text = element_text(size = legend_size),
            legend.position = legend_position)
}

#' plot_stacked_bar
#'
#' This is a function for stacked barplot.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#' metadata and phylogenetic tree. Set it to NULL if using an OTU table and
#' metadata to draw the stacked barplot.
#' @param feature The feature that shows in x-axis text with different colors.
#' Default is NA, x-axis will show in black.
#' @param level  A taxonomy level to plot stacked bar. Default is NA. Required
#' when using phyloseq object.
#' @param x_size The front size of x-axis text. Default is 8.
#' @param legend_position Legend position. Default is "top". legend_position
#' should be one of "none", "left", "right", "bottom", "top".
#' @param legend_size Lengend size. Default is 10.
#' @param order A vector of arranged SampleID in wanted order. Default is NULL.
#' If NULL, plot_stacked_bar will automatically arrange samples by the most
#' abundance taxonomy in decreasing order.
#' @param otu_table An OTU table. Taxa are colnames, SampleID are rownames
#' (just like phyloseq). Set it to NULL if using phyloseq.
#' @param metadata A metadata for the OTU table. SampleID are rownames (just
#' like phyloseq). Set it to NULL if using phyloseq.
#' @param relative_abundance Turn plot into relative abundance or not. Default
#' is FALSE.
#' @export
#' @examples
#' plot_stacked_bar(demo_phyloseq_object, level = "Order")

plot_stacked_bar <- function(phyloseq = NULL, level = NA, feature = NA,
                             otu_table = NULL, metadata = NULL, order = NULL,
                             relative_abundance = FALSE, x_size = 8,
                             legend_position = "top", legend_size = 10) {
  # Detact variables
  if (is.null(phyloseq)) {
    if (is.null(otu_table)) {
      # If no phyloseq, require otu_table
      stop("Argument 'phyloseq' and 'otu_table' are both missing,
           please input a phyloseq object or an OTU table and metadata.")
    } else {
      if (is.null(metadata)) {
        # If have otu_table, require metadata
        stop("Argument 'metadata' is missing, please input a metadata for the
             OTU table.")
      }
    }
  } else {
    if (!is.null(otu_table)) {
      # If have phyloseq, require no otu_table
      stop("Argument 'phyloseq' and 'otu_table' are both detected, please input
           one of them, not both.")
    } else if (is.na(level)) {
      # If have phyloseq, require level
      stop("Argument 'level' is missing. Plaese choose a taxonomy level to plot
           stacked bar.")
    }
  }
  # First construct otu then convert to percentage
  if (!is.null(phyloseq)) {
    otu <- construct_otu_table(phyloseq, level)
  } else {
    otu <- otu_table
  }
  if (relative_abundance) {
    otu <- otu %>% convert_to_percentage(row_sum = TRUE)
  }
  # Construct table for stacked bar plot
  plot_tab <- otu %>%
    # First re-order taxonomy by total counts
    .[,order(colSums(.), decreasing = TRUE)] %>%
    # Then re-order sample by the most abundance taxonomy
    .[order(.[,1], decreasing = TRUE),] %>%
    # Turn SampleID to a column
    rownames_to_column(var = "SampleID")
  # Prepare levels for taxonomy levels
  levels_level <- colnames(plot_tab)[2:ncol(plot_tab)]
  # Extract metadata
  if (is.null(otu_table)) {
    sample_feature <- extract_metadata_phyloseq(phyloseq)
  } else {
    sample_feature <- metadata %>% rownames_to_column("SampleID")
  }
  # Prepare levels for SampleID
  if (is.null(order)) {
    levels_SampleID <- plot_tab$SampleID
  } else if (length(order) != length(plot_tab$SampleID)) {
    stop("The length of order and SampleID does not match.")
  } else if (all(order %in% plot_tab$SampleID)) {
    if (all(plot_tab$SampleID %in% order)) {
      levels_SampleID <- order
    }
  } else {
    stop("Given order must contain all SampleID.")
  }
  # Join metadata
  plot_tab <- left_join(plot_tab, sample_feature)
  # Add levels to SampleID
  plot_tab$SampleID <- factor(plot_tab$SampleID, levels = levels_SampleID)
  # Prepare levels for feature (colors in x-axis)
  if (!is.na(feature)) {
    if (is.null(order)) {
      levels_feature <- plot_tab[[feature]] %>% as.factor()
    } else {
      levels_feature <- as.data.frame(order)
      colnames(levels_feature) <- "Sort"
      levels_feature <- levels_feature %>%
        # Arrange feature order by 'order' variable
        left_join(select(plot_tab, SampleID, !!feature),
                  by = c("Sort" = "SampleID")) %>%
        # Drop levels
        .[[2]] %>% as.character() %>% as.factor()
    }
  } else {
    levels_feature <- "black"
  }
  # Turn plot_table to a long table for plotting
  plot_tab <- gather(plot_tab,
                     colnames(plot_tab)[2:(ncol(plot_tab)-1)],
                     key = level,
                     value = "abundance")
  # Add levels to taxonomy levels
  plot_tab$level <- factor(plot_tab$level, levels = levels_level)
  # Bar plot
  ggplot(plot_tab, aes(x = SampleID, y = abundance)) +
    scale_fill_manual(values = distinctive_colors) +
    geom_bar(mapping = aes(fill = level), stat = "identity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.text.x = element_text(size = x_size, angle = 90,
                                     # 'color' can only work with factor
                                     color = levels_feature),
          legend.text = element_text(size = legend_size),
          legend.position = legend_position)
}

#' plot_stacked_bar
#'
#' This is a function for stacked bar plot.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#'                 metadata and phylogenetic tree. If you don't have a phyloseq
#'                 object, you can input OTU table and metadata to draw stacked
#'                 barplot.
#' @param feature The feature that shows in x-axis text with different colors.
#' @param level  A taxonomy level to plot stacked bar.
#' @param x_size The front size of x-axis text. Default is 8.
#' @param legend_position Legend position. Default is "top". legend_position
#'                        should be one of "none", "left", "right", "bottom",
#'                        "top".
#' @param legend_size Lengend size. Default is 10.
#' @param order A vector of arranged SampleID in wanted order. Default is NA.
#'              If NA, plot_stacked_bar will automatically arrange samples by
#'              the most abundance taxonomy in decreasing order.
#' @param otu_table An OTU table. Taxa are colnames, SampleID are rownames
#'                  (just like phyloseq).
#' @param metadata A metadata for the OTU table. SampleID are rownames (just
#'                 like phyloseq).
#' @export
#' @examples
#' plot_stacked_bar(demo_phyloseq_object, feature = "diagnosis",
#'                  level = "Order")

plot_stacked_bar <- function(phyloseq = NA, level = NA, feature = NA,
                             order = NA, x_size = 8, legend_position = "top",
                             legend_size = 10, otu_table = NA, metadata = NA) {
  # Detact variables
  if (is.na(phyloseq)) {
    if (is.na(otu_table)) {
      # If no phyloseq, require otu_table
      stop("Argument 'phyloseq' and 'otu_table' are both missing,
           please input a phyloseq object or an OTU table and metadata.")
    } else {
      if (is.na(metadata)) {
        # If have otu_table, require metadata
        stop("Argument 'metadata' is missing, please input a metadata for the
             OTU table.")
      }
    }
  } else {
    if (!is.na(otu_table)) {
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
  if (!is.na(phyloseq)) {
    otu <- construct_otu_table(phyloseq, level)
  } else {
    otu <- otu_table
  }
  otu_percent <- otu %>%
    convert_to_percentage(row_sum = TRUE) %>%
    as.data.frame()
  # Construct table for stacked bar plot
  plot_tab <- otu_percent %>%
    # First re-order taxonomy by total counts
    .[,order(colSums(.), decreasing = TRUE)] %>%
    # Then re-order sample by the most abundance taxonomy
    .[order(.[,1], decreasing = TRUE),] %>%
    # Turn SampleID to a column
    rownames_to_column(var = "SampleID")
  # Prepare levels for taxonomy levels
  levels_level <- colnames(plot_tab)[2:ncol(plot_tab)]
  # Extract metadata
  if (is.na(otu_table)) {
    sample_feature <- extract_metadata_phyloseq(phyloseq, feature)
  } else {
    sample_feature <- metadata %>% rownames_to_column("SampleID")
  }
  # Prepare levels for SampleID
  if (all(is.na(order))) {
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
    if (all(is.na(order))) {
      levels_feature <- plot_tab[,ncol(plot_tab)] %>% as.factor()
    } else {
      levels_feature <- as.data.frame(order) %>%
        # Arrange feature order by 'order' variable
        left_join(plot_tab[,c(1, ncol(plot_tab))],
                  by = c("order" = "SampleID")) %>%
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
    geom_bar(mapping = aes(fill = level), position = "fill",
             stat = "identity") +
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

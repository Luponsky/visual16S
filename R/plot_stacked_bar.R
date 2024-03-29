#' plot_stacked_bar
#'
#' This is a function for stacked barplot.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample metadata and phylogenetic
#' tree. Set it to NULL if using an OTU table and metadata to draw the stacked barplot.
#'
#' @param feature The feature that will show as a colorbar under x axis. Default is NA.
#'
#' @param level A taxonomy level to plot stacked bar. Default is NA. Required when using phyloseq
#' object.
#'
#' @param legend_position Legend position. Default is "top". legend_position should be one of "none",
#' "left", "right", "bottom", "top".
#'
#' @param order A vector of arranged SampleID in wanted order. Default is NULL. If NULL,
#' plot_stacked_bar will automatically arrange samples by the most abundance taxonomy in decreasing
#' order.
#'
#' @param otu_table An OTU table. Taxa are colnames, SampleID are rownames (just like phyloseq). Set it
#' to NULL if using phyloseq.
#'
#' @param metadata A metadata for the OTU table. SampleID are rownames (just like phyloseq). Set it to
#' NULL if using phyloseq.
#'
#' @param relative_abundance Turn plot into relative abundance or not. Default is FALSE.
#'
#' @param colors A vector of colors. The number of colors should be larger than the number of level.
#' Default is NULL, if NULL, plot_stacked_bar will use distinctive_colors for the plot.
#'
#' @export
#'
#' @examples
#' plot_stacked_bar(demo_phyloseq_object, level = "Order")

plot_stacked_bar <- function (
  phyloseq = NULL,
  level = NA,
  feature = NA,
  otu_table = NULL,
  metadata = NULL,
  order = NULL,
  relative_abundance = FALSE,
  legend_position = "top",
  colors = NULL
) {
  # Detact variables
  if (is.null(phyloseq)) {
    if (is.null(otu_table)) {
      # If no phyloseq, require otu_table
      stop(paste0("Argument 'phyloseq' and 'otu_table' are both missing, please input a phyloseq ",
                  "object or an OTU table with metadata."))
    } else {
      if (is.null(metadata)) {
        # If have otu_table, require metadata
        stop(paste0("Argument 'metadata' is missing, please input a metadata for the OTU table."))
      }
    }
  } else {
    if (!is.null(otu_table)) {
      # If have phyloseq, require no otu_table
      stop(paste0("Argument 'phyloseq' and 'otu_table' are both detected, please input one of them, ",
                  "not both."))
    } else if (is.na(level)) {
      # If have phyloseq, require level
      stop(paste0("Argument 'level' is missing. Plaese choose a taxonomy level to plot stacked bar."))
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
    # Construct table for stacked bar plot
    plot_tab <- otu %>%
      # First re-order taxonomy by total counts
      .[,order(colSums(.), decreasing = TRUE)] %>%
      # Then re-order sample by the most abundance taxonomy
      .[order(.[,1], decreasing = TRUE),] %>%
      # Turn SampleID to a column
      rownames_to_column(var = "SampleID")
    y_label <- "Relative Abundance"
  } else {
    plot_tab <- otu %>%
      # First re-order taxonomy by total counts
      .[,order(colSums(.), decreasing = TRUE)] %>%
      # Then re-order sample by total abundance
      .[order(rowSums(.), decreasing = TRUE),] %>%
      # Turn SampleID to a column
      rownames_to_column(var = "SampleID")
    y_label <- "Raw Count"
  }
  # Prepare levels for taxonomy levels (order within each bar)
  levels_level <- colnames(plot_tab)[2:ncol(plot_tab)]
  # Extract metadata
  if (is.null(otu_table)) {
    sample_feature <- extract_metadata_phyloseq(phyloseq)
  } else {
    sample_feature <- metadata %>% rownames_to_column("SampleID")
  }
  # Prepare levels for SampleID (order of x.axis)
  if (is.null(order)) {
    #if (!is.na(feature)) {
    #  # Save origin colnames
    #  original_colnames <- colnames(plot_tab)
    #  # Add a new column that is rowsum
    #  plot_tab$total <- plot_tab %>%
    #    column_to_rownames("SampleID") %>%
    #    rowSums()
    #  plot_tab <- plot_tab %>%
    #    # Join metadata
    #    left_join(sample_feature) %>%
    #    # Group by feature parameter
    #    group_by_(feature) %>% # group_by_() can pass variable to group_by()
    #    # Arrange rows by total abundance and by group
    #    arrange(desc(total), .by_group = TRUE) %>% # Set .by_group = TRUE to arrange by group
    #    select(original_colnames)
    #  levels_SampleID <- plot_tab$SampleID
    #} else {
      levels_SampleID <- plot_tab$SampleID
    #}
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
  ## Prepare levels for feature (colors of x.axis)
  ## 要想X轴或Y轴显示颜色，必须添加一组颜色的vector，顺序和X轴或Y轴顺序一致，vector内容是不同的level
  ## Add "theme(axis.text.x = element_text(color = levels_feature))" to ggplot.
  ## "color" parameter can only work with factor.
  #if (!is.na(feature)) {
  #  if (is.null(order)) {
  #    levels_feature <- plot_tab[[feature]] %>% as.factor()
  #  } else {
  #    levels_feature <- as.data.frame(order)
  #    colnames(levels_feature) <- "Sort"
  #    levels_feature <- levels_feature %>%
  #      # Arrange feature order by 'order' variable
  #      left_join(select(plot_tab, SampleID, !!feature), by = c("Sort" = "SampleID")) %>%
  #      # Drop levels
  #      .[[2]] %>% as.character() %>% as.factor()
  #  }
  #} else {
  #  levels_feature <- "black"
  #}
  # Turn plot_table to a long table for plotting
  plot_tab <- gather(plot_tab, levels_level, key = level, value = "abundance")
  # Add levels to taxonomy levels
  plot_tab$level <- factor(plot_tab$level, levels = levels_level)
  # Bar plot
  p1 <- ggplot(plot_tab, aes(x = SampleID, y = abundance)) +
    geom_bar(mapping = aes(fill = level), stat = "identity") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 14, angle = 90),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.position = legend_position,
      legend.title = element_blank()
    ) +
    ylab(y_label)
  if (is.null(colors)) {
    p1 <- p1 + scale_fill_manual(values = distinctive_colors)
  } else if (length(colors) < length(levels_level)) {
    stop(paste0("The number of colors is smaller than the number of ", level, "."))
  } else {
    p1 <- p1 + scale_fill_manual(values = colors)
  }
  if (is.na(feature)) {
    p1
  } else {
    # If feature is given, make abundance plot (p1) & feature colorbar (p2) and align them
    p1 <- p1 +
      theme(
        #panel.grid = element_blank(), # Remove grid line under the plot
        #panel.border = element_blank(), # Remove border
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = legend_position,
        legend.title = element_blank()
      )
    x = "SampleID"
    y = 1
    fill = feature
    p2 <- ggplot(data = plot_tab, aes_string(x = x, y = y, fill = fill)) +
      geom_tile() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 14, angle = 90),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom"
      ) +
      scale_fill_manual(values = distinctive_colors[4:length(distinctive_colors)])
    ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(8, 1.5), align = "v")
  }
}

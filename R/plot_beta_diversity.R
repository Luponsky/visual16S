#' plot_beta_diversity
#'
#' This is a function for plotting beta diversity.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample metadata and phylogenetic
#' tree.
#'
#' @param feature The column name of the feature you want to select from metadata, e.g. "Phenotype".
#'
#' @param feature2 The column name of another feature you want to select from metadata, which will show
#' in different shape, e.g. "Gender". Default is NA.
#'
#' @param method The method to calculate beta diversity. Method should be one of c("bray", "jaccard",
#' "unifrac", "wunifrac"). Default is "bray". PS: "unifrac" and "wunifrac" require a phylogenetic tree.
#'
#' @param colors A vector of colors. The number of colors should be larger than the number of feature.
#' Default is NULL, if NULL, plot_beta_diversity will use ggsci::scale_color_jco for the plot.
#'
#' @export
#'
#' @examples
#' plot_beta_diversity(demo_phyloseq_object, feature = "diagnosis")

plot_beta_diversity <- function(
  phyloseq,
  feature,
  feature2 = NA,
  method = "bray",
  colors = NULL
  ){
  set.seed(99)
  ## Step 1: Calculate beta diversity
  if (!method %in% c("bray", "jaccard", "unifrac", "wunifrac")) {
    stop(paste0('Beta diversity method should be one of c("bray", "jaccard", "unifrac", "wunifrac").'))
  } else if (method %in% c("unifrac", "wunifrac")) {
    # Requires phyloseq-class that contains both an otu_table and a phylogenetic tree
    beta <- cmdscale(phyloseq::distance(physeq = phyloseq, method = method), eig = TRUE)
  } else {
    beta <- cmdscale(vegan::vegdist(otu_table(phyloseq), method = method), eig = TRUE)
  }
  ## Step 2: Construct plot table
  # Extract PC1 and PC2
  PC <- as.data.frame(beta$points)
  colnames(PC) <- c("PC1", "PC2")
  PC <- rownames_to_column(PC, var = "SampleID")
  # Extract metadata
  metadata <- extract_metadata_phyloseq(phyloseq)
  # Join two tables
  beta_plot <- left_join(PC, metadata)
  # Print beta-diversity table
  if (is.na(feature2)) {
    select(beta_plot, SampleID, !!feature, PC1, PC2) %>% print()
  } else {
    select(beta_plot, SampleID, !!feature, !!feature2, PC1, PC2) %>% print()
  }
  ## Step 3: Plot beta diversity
  # Make x-axis and y-axis names for aes_string
  x_name <- "PC1"
  y_name <- "PC2"
  # Make factor for color, prevent "Error: Continuous value supplied to discrete scale"
  beta_plot[[feature]] <- beta_plot[[feature]] %>% as.factor()
  # Plot beta diversity
  if (is.na(feature2)) {
    p <- ggplot(data = beta_plot,
                # aes_string() can pass variables to ggplot, aes() can't
                aes_string(x = x_name, y = y_name, color = feature)) +
      geom_point(size = 3) +
      xlab(paste("PC1:", round(100*as.numeric(beta$eig[1]/sum(beta$eig)), 2), "%", sep = " ")) +
      ylab(paste("PC2:", round(100*as.numeric(beta$eig[2]/sum(beta$eig)), 2), "%", sep = " ")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 8))
  } else {
    p <- ggplot(data = beta_plot,
                # Use aes_string() to pass variables to ggplot
                aes_string(x = x_name, y = y_name, color = feature, shape = feature2)) +
      geom_point(size = 3) +
      xlab(paste("PC1:", round(100*as.numeric(beta$eig[1]/sum(beta$eig)), 2), "%", sep = " ")) +
      ylab(paste("PC2:", round(100*as.numeric(beta$eig[2]/sum(beta$eig)), 2), "%", sep = " ")) +
      scale_shape_manual(values = c(0:6)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 8))
  }
  if (is.null(colors)) {
    p + ggsci::scale_color_jco() + ggsci::scale_fill_jco()
  } else if (length(colors) < length(unique(beta_plot[[feature]]))) {
    stop(paste0("The number of colors is smaller than the number of ", feature, "."))
  } else {
    p + scale_color_manual(values = colors)
  }
}

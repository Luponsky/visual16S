#' plot_beta_diversity
#'
#' This is a function for plotting beta diversity.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#'                 metadata and phylogenetic tree.
#' @param feature The column name of the feature you want to select from
#'                metadata, e.g. "Phenotype".
#' @param feature2 The column name of another feature you want to select from
#'                 metadata, which will show in different shape, e.g. "Gender".
#'                 Default is NA.
#' @param method The method to calculate beta diversity. Method should be one
#'               of "bray", "jaccard", "unifrac", "wunifrac".
#' @export
#' @examples
#' plot_beta_diversity(demo_phyloseq_object, feature = "diagnosis",
#'                     method = "bray")

plot_beta_diversity <- function(phyloseq, feature, feature2 = NA, method){
  set.seed(99)
  ## Step 1: Calculate beta diversity
  if (!method %in% c("bray", "jaccard", "unifrac", "wunifrac")) {
    stop('beta diversity method should be one of "bray", "jaccard", "unifrac", "wunifrac".')
  } else if (method %in% c("unifrac", "wunifrac")) {
    # Requires phyloseq-class that contains both an otu_table and a
    # phylogenetic tree
    beta <- cmdscale(phyloseq::distance(physeq = phyloseq, method = method),
                     eig = TRUE)
  } else {
    beta <- cmdscale(vegan::vegdist(otu_table(phyloseq), method = method),
                     eig = TRUE)
  }
  ## Step 2: Construct plot table
  # Extract PC1 and PC2
  PC <- as.data.frame(beta$points)
  colnames(PC) <- c("PC1", "PC2")
  PC <- rownames_to_column(PC, var = "SampleID")
  # Print beta-diversity table
  print(PC)
  # Extract metadata
  metadata <- extract_metadata_phyloseq(phyloseq)
  # Join two tables
  beta_plot <- left_join(PC, metadata)
  ## Step 3: Plot beta diversity
  # Make x-axis and y-axis names for aes_string
  x_name <- "PC1"
  y_name <- "PC2"
  # Make factor for color, prevent "Error: Continuous value supplied to
  # discrete scale"
  beta_plot[[feature]] <- beta_plot[[feature]] %>% as.factor()
  # Plot beta diversity
  if (is.na(feature2)) {
    p <- ggplot(data = beta_plot,
                # aes_string() can pass variables to ggplot, aes() can't
                aes_string(x = x_name, y = y_name, color = feature)) +
      geom_point(size = 3) +
      xlab(paste("PC1:", round(100*as.numeric(beta$eig[1]/sum(beta$eig)), 2),
                 "%", sep = " ")) +
      ylab(paste("PC2:", round(100*as.numeric(beta$eig[2]/sum(beta$eig)), 2),
                 "%", sep = " ")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 14))
    p + ggsci::scale_color_jco() + ggsci::scale_fill_jco()
  } else {
    p <- ggplot(data = beta_plot,
                # Use aes_string() to pass variables to ggplot
                aes_string(x = x_name, y = y_name, color = feature,
                           shape = feature2)) +
      geom_point(size = 3) +
      xlab(paste("PC1:", round(100*as.numeric(beta$eig[1]/sum(beta$eig)), 2),
                 "%", sep = " ")) +
      ylab(paste("PC2:", round(100*as.numeric(beta$eig[2]/sum(beta$eig)), 2),
                 "%", sep = " ")) +
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

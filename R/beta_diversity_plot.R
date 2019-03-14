#' beta_diversity_plot
#'
#' This is a function for plotting beta diversity.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample metadata and
#'                 phylogenetic tree.
#' @param feature The column name of the feature you want to select from metadata, e.g. "Phenotype".
#' @param feature2 The column name of another feature you want to select from metadata, e.g. "Gender".
#'                 Default is NA.
#' @param method The method to calculate beta diversity. Method should be one of "bray", "jaccard",
#'               "unifrac", "wunifrac".
#' @param size The size of the plot. Default is 3.
#' @export
#' @examples
#' beta_diversity_plot(Shaoyifu_phyloseq, feature = "diagnosis", method = "bray")

beta_diversity_plot <- function(phyloseq, feature, feature2 = NA, method, size = 3){
  ## Step 1: Calculate beta diversity
  if (!method %in% c("bray", "jaccard", "unifrac", "wunifrac")) {
    stop('beta diversity method should be one of "bray", "jaccard", "unifrac", "wunifrac".')
  } else if (method %in% c("unifrac", "wunifrac")) {
    beta_diversity <- cmdscale(phyloseq::distance(physeq = phyloseq, method = method), eig = TRUE)
  } else {
    beta_diversity <- cmdscale(vegan::vegdist(otu_table(phyloseq), method = method), eig = TRUE)
  }
  ## Step 2: Construct plot table
  # Extract PC1 and PC2
  PC <- as.data.frame(beta_diversity$points)
  colnames(PC) <- c("PC1", "PC2")
  PC <- rownames_to_column(PC, var = "subject_id")
  # Extract metadata
  metadata <- extract_metadata_from_phyloseq(phyloseq)
  # Mrtge
  beta_diversity_tab <- left_join(PC, metadata)
  ## Step 3: Plot beta diversity
  # Make x-axis and y-axis names for aes_string
  x_name <- "PC1"
  y_name <- "PC2"
  # Make factor for color, prevent "Error: Continuous value supplied to discrete scale"
  beta_diversity_tab[[feature]] <- beta_diversity_tab[[feature]] %>% as.factor()
  # Plot beta diversity
  # aes_string() can pass variables to ggplot, aes() can't
  if (is.na(feature2)) {
    ggplot(data = beta_diversity_tab, aes_string(x = x_name, y = y_name, color = feature)) +
      geom_point(size = size) +
      xlab(paste("PC1:",
                 round(100*as.numeric(beta_diversity$eig[1]/sum(beta_diversity$eig)), 2),
                 "%",
                 sep = " ")) +
      ylab(paste("PC2:",
                 round(100*as.numeric(beta_diversity$eig[2]/sum(beta_diversity$eig)), 2),
                 "%",
                 sep = " ")) +
      scale_color_manual(values = c("#00AFBB", "#FC4E07", "#7FC97F", "#BEAED4"))
  } else {
    ggplot(data = beta_diversity_tab, aes_string(x = x_name, y = y_name, color = feature,
                                                 shape = feature2)) +
      geom_point(size = size) +
      xlab(paste("PC1:",
                 round(100*as.numeric(beta_diversity$eig[1]/sum(beta_diversity$eig)), 2),
                 "%",
                 sep = " ")) +
      ylab(paste("PC2:",
                 round(100*as.numeric(beta_diversity$eig[2]/sum(beta_diversity$eig)), 2),
                 "%",
                 sep = " ")) +
      scale_color_manual(values = c("#00AFBB", "#FC4E07", "#7FC97F", "#BEAED4")) +
      scale_shape_manual(values = c(0:6))
  }
}

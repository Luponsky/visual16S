#' plot_sparsity
#'
#' plot_sparsity can show the sparsity of an OTU table. It will calculate the
#' prevalence of each OTU in all samples. X-axis is an OTU's prevalence. Y-axis
#' is OTU count.
#'
#' @param otu_table An OTU table in phyloseq object format. Rownames are
#' sample ID, colnames are taxa.
#'
#' @param binwidth The width of the bins. Can be specified as a numeric value,
#' or a function that calculates width from x. The default is to use bins bins
#' that cover the range of the data. You should always override this value,
#' exploring multiple widths to find the best to illustrate the stories in your
#' data.
#'
#' @export
#'
#' @examples
#' plot_sparsity(dada2_res$seq_tab, 5)

plot_sparsity <- function(otu_table, binwidth = NA) {
  otu_table <- otu_table %>%
    # Protect rownames (Tidyverse will automatically remove rownames)
    rownames_to_column() %>%
    # Remove all 0 OTUs
    filter_if(is.numeric, any_vars(. != 0)) %>%
    # Recover rownames
    column_to_rownames() %>%
    # Transposes OTU table
    t() %>% as.data.frame()
  # Replace 0 with NA
  otu_table[otu_table == 0] <- NA
  # Calculate sparsity
  otu_table <- apply(
    otu_table, 1, function(x) round((sum(!is.na(x))/ncol(otu_table))*100, 0)
    ) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'otu')
  colnames(otu_table)[2] <- 'prevalence'
  if (!is.na(binwidth)) {
    ggplot(otu_table, aes(prevalence)) +
      geom_histogram(binwidth = binwidth) +
      xlab('Prevalence of each OTU') +
      ylab('Count') +
      theme_bw() +
      theme(panel.grid = element_blank())
  } else {
    ggplot(otu_table, aes(prevalence)) +
      geom_histogram() +
      xlab('Prevalence of each OTU') +
      ylab('Count') +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
}

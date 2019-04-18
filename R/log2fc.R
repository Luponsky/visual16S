#' log2fc
#'
#' This is a function for plotting log2 fold change.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample
#'                 metadata and phylogenetic tree.
#' @param feature The column name of the feature you want to select from metadata,
#'                e.g. "Phenotype".
#' @param level Which taxonomy level to calculate fold change. Default is NA. If
#'              level is given, will use construct_otu_table function to construct
#'              OTU table, and use DESeq to calculate fold change.
#'
#' @param p_value The cut off P value for the fold change. Default is 0.01.
#' @param save_res Default is FALSE. If TRUE, will save original result
#'                 DESeq2_result.rds to current working directory.
#' @param reference The control group. Default is NA.
#' @param treatment The treatment group. Default is NA.
#' @export
#' @examples
#' log2fc(demo_phyloseq_object, feature = "diagnosis", level = "Genus")

log2fc <- function(phyloseq, feature, level = NA, p_value = 0.01,
                   save_res = FALSE, reference = NA, treatment = NA) {
  set.seed(99)
  ## Step 1: Construt table for DESeq2
  # Create a string to parse feature argument to DESeq
  feature_formula <- paste0("~ ", feature)
  # Create a DESeq Data Set
  if (is.na(level)) {
    dds <- phyloseq_to_deseq2(phyloseq, design = as.formula(feature_formula))
    cts <- counts(dds)
    geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0
                      else exp(sum(log(row[row != 0]))/length(row)))
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  } else {
    countData <- construct_otu_table(phyloseq, level) %>% t()
    colData <- extract_metadata_phyloseq(phyloseq = phyloseq, feature = feature)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData,
                                  design = as.formula(feature_formula))
  }
  dds <- DESeq(dds)
  ## Step 2: Calculate log2 fold change
  if (!is.na(reference)) {
    res <- results(dds, contrast = c(feature, treatment, reference))
  } else {
    res <- results(dds)
  }
  res <- res[order(res$padj),]
  ## Step 3: Plot fold change
  # Create table for plotting
  ifelse(is.na(level), var <- "OTU", var <- level)
  log2fc <- res %>% as.data.frame() %>% rownames_to_column(var = var) %>%
    filter(!is.na(padj))
  # Identify differential expression
  log2fc$significant = ifelse(log2fc$padj > p_value,
                              paste0('p-value > ', p_value),
                              paste0('p-value < ', p_value))
  if (save_res) {
    # save res object to current working directory
    saveRDS(log2fc, "DESeq2_result.rds")
    print('DESeq2_result.rds has been saved to current working directory.')
  }
  # Extract table that meet the p-value threshold to print to the screen
  log2fc_print <- log2fc %>%
    .[.$padj < p_value,] %>%
    select(!!var, log2FoldChange, padj) %>%
    .[order(.[,2], decreasing = TRUE),]
  # EnhancedVolcano plot
  print(res@elementMetadata$description[2])
  print(log2fc_print)
  EnhancedVolcano(log2fc, lab = log2fc[[var]], x = "log2FoldChange", y = "padj",
                  FCcutoff = 1, pCutoff = p_value,
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  transcriptPointSize = 3,
                  col = c("darkgrey", "#00AFBB", "#FC4E07", "red2"),
                  legend = c("Not Significant",
                             paste0("p > ", p_value, ", log2 fold change > 1"),
                             paste0("p < ", p_value, ", log2 fold change < 1"),
                             paste0("p < ", p_value, ", log2 fold change > 1"))
  )
}

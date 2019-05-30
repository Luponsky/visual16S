#' construct_lefse_table
#'
#' construct_lefse_table can construct a LEfSe-format OTU table.
#' See http://huttenhower.sph.harvard.edu/galaxy/ -> LEfSe -> Format Data for LEfSe for more details.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, samplemetadata and phylogenetic
#' tree.
#'
#' @param feature The column name of the feature you want to select. In final table, feature will be the
#' first row.
#'
#' @param level The coloumn name of the level wanted to select. Default is "all". If "all" then retain
#' all taxonomy level, else retain the taxonomy from Kingdom to selected level, drop everything else.
#' Level name should be one of c("all", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
#' "Species"). Taxonomy will be seperated by "|".
#'
#' @export
#'
#' @examples
#' construct_lefse_table(demo_phyloseq_object, feature = "diagnosis", level = "all") %>% .[1:20,1:5]

construct_lefse_table <- function(phyloseq, feature, level = "all") {
  # Set options, prevent R turnning numeric value to factor
  options(stringsAsFactors = FALSE)
  # Read in phyloseq object
  # Convert OTU to relative abundance
  otu <- otu_table(phyloseq) %>% as.data.frame() %>%
    convert_to_percentage() %>%
    t() %>% as.data.frame() %>%
    rownames_to_column(var = "OTU_ID")
  taxa <- tax_table(phyloseq) %>% as.data.frame() %>%
    rownames_to_column(var = "OTU_ID")
  # Select taxa by given level
  if (level == "all") {
    levels <- colnames(taxa)[2:ncol(taxa)]
  } else {
    taxa <- taxa %>% select(OTU_ID:!!level)
    levels <- colnames(taxa)[2:ncol(taxa)]
  }
  # Extract levels name
  levels <- colnames(taxa)[2:ncol(taxa)]
  # Remove rows in taxa that is all NA
  all_NA_taxa <- taxa %>% filter_all(all_vars(is.na(.)))
  taxa <- taxa %>% filter(!OTU_ID %in% all_NA_taxa$OTU_ID)
  # Convert metadata and create lefse table
  metadata <- extract_metadata_phyloseq(phyloseq = phyloseq, feature = feature)
  lefse <- metadata %>% t() %>% as.data.frame()
  colnames(lefse) <- lefse %>% .[rownames(.) == "SampleID",]
  # Combine otu to lefse table
  for(i in levels) {
    taxa_tmp <- taxa %>%
      as.data.frame() %>%
      select(1:which(colnames(taxa) == i)) %>%
      filter_all(all_vars(!is.na(.)))
    taxa_tmp <- taxa_tmp %>% unite(Taxonomy, 2:ncol(taxa_tmp), sep = "|")
    otu_tmp <- otu %>% filter(OTU_ID %in% taxa_tmp$OTU_ID) %>%
      left_join(taxa_tmp)
    otu_tmp <- otu_tmp %>%
      #group_by_() can pass variable to goup_by() function
      group_by(Taxonomy) %>%
      summarise_if(is.numeric, sum, na.rm=TRUE) %>%
      column_to_rownames(var = "Taxonomy")
    lefse <- rbind(lefse, otu_tmp)
  }
  lefse <- rownames_to_column(lefse)
  colnames(lefse) <- 1:ncol(lefse)
  return(lefse)
}

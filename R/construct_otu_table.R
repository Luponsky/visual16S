#' construct_otu_table
#'
#' construct_otu_table can construct a OTU table with a phyloseq object.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample metadata and
#'                 phylogenetic tree.
#' @param level The coloumn name of the level wanted to select. Default is "all". If "all" then retain
#'              all taxonomy level and seperate by "; ", else ONLY retain the given taxonomy level,
#'              drop everything else. Level name should be one of "all", "Kingdom", "Phylum", "Class",
#'              "Order", "Family", "Genus", "Species".
#' @export
#' @examples
#' construct_otu_table(Shaoyifu_phyloseq, level = "Order") %>% .[,1:5]

construct_otu_table <- function(phyloseq, level = "all") {
  # Set options, prevent R turnning numeric value to factor
  options(stringsAsFactors = FALSE)
  # Check if input 'level' is correct
  if (!level %in% c("all", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    stop('level should be one of "all", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species".')}
  # Read in sequence table and taxonomy table from phyloseq
  otu <- otu_table(phyloseq) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    rownames_to_column(var = "OTU_ID")
  taxa <- tax_table(phyloseq) %>% as.data.frame() %>%
    rownames_to_column(var = "OTU_ID")
  ## Step 1: Clean sequence table and taxonomy table
  # If a sequence's taxonomy is all NA, remove it
  all_NA_taxa <- taxa %>% filter_all(all_vars(is.na(.)))
  taxa <- taxa %>% filter(!OTU_ID %in% all_NA_taxa$OTU_ID)
  ## Step 2: Construct OTU table
  if (level == "all") {
    # If level == "all", collapse all taxonomy levels and separate by "; "
    level <- "Taxonomy"
    taxa <- taxa %>% unite(Taxonomy, 2:ncol(taxa), sep = "; ")
  } else {
    # If level != "all", select the given taxonomy level
    taxa <- taxa %>% select(OTU_ID, level) %>% filter(., !is.na(.[[level]]))
    otu <- otu %>% filter(OTU_ID %in% taxa$OTU_ID)
  }
  ## Step 3: Merge sequence table and taxonomy table
  otu <- left_join(otu, taxa)
  ## Step 4: Add abundance of those have the same taxonomy names and convert it to data frame
  otu <- otu %>%
    group_by_(level) %>% #group_by_() can pass variable to goup_by() function
    summarise_if(is.numeric, sum, na.rm=TRUE) %>%
    filter_if(is.numeric, any_vars(. != 0)) %>%
    column_to_rownames(var = level) %>%
    t() %>%
    as.data.frame()
  return(otu)
}

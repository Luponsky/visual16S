#' extract_metadata_from_phyloseq
#'
#' extract_metadata_from_phyloseq can extract metadata from a phyloseq object. First, the function
#' will extract metadata from phyloseq object using phyloseq::sample_data, and turn it into a tibble
#' which will turn rownames into a column name 'subject_id'. If feature parameter is given, then will
#' select subject_id and feature column, and add levels to the selected feature column.
#'
#' @param phyloseq A phyloseq object contain otu table, taxonomy table, sample metadata and
#'                 phylogenetic tree.
#' @param feature The column name of the feature you want to select. Default is NA. If NA, will return
#'                the complete metadata, else will return subject id and feature column that's given.
#' @export
#' @examples
#' extract_metadata_from_phyloseq(Shaoyifu_phyloseq, feature = "diagnosis")

extract_metadata_from_phyloseq <- function(phyloseq, feature = NA) {
  ## Step 1: Extract metadata from phyloseq and turn into tibble
  metadata <- phyloseq %>%
    sample_data() %>%
    as.matrix() %>%
    as.data.frame()
  if ("subject_id" %in% colnames(metadata)) {
    select(metadata, -subject_id) %>%
      rownames_to_column(var = "subject_id")
    warning("replace 'subject_id' column with rownames")
  } else {
    metadata <- rownames_to_column(metadata, var = "subject_id")
  }
  if (is.na(feature)) {
    return(metadata)
  } else {
    # Select column by feature name
    metadata <- dplyr::select(metadata, subject_id, !!feature)
    # Add levels to column
    metadata[[feature]] <- metadata[[feature]] %>% as.character() %>% factor()
    return(metadata)
  }
}

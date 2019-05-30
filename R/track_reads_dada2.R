#' track_reads_dada2
#'
#' track_reads_dada2 can track the reads count in a dada2 workflow result which created by Xbiome 16S
#' pipeline. Xbiome 16S pipeline dada2 workflow will generate a list that contain sequence table,
#' taxonomy table and reads track data frame. Input the reads track data frame and read type, this
#' function can draw a line plot of reads track of every sample. X-axis will be every stage in dada2
#' workflow, Y-axis will be the reads counts.
#'
#' @param reads_track The reads track data frame from Xbiome 16S pipeline dada2 workflow result.
#'
#' @param single_end Default is FALSE. If single_end == TRUE, means the sequence files are single end,
#' the x-axis will contain 'input', 'filtered', 'dereplicated', 'nonchim'. If single_end == FALSE,
#' means the sequence files are paired end, the x-axis will contain 'input', 'filtered', 'denoisedF',
#' 'denoisedR', 'merged', 'nonchim'.
#'
#' @param relative_abundance Default is FALSE. If TRUE, will turn values to relative abundance.
#'
#' @param legend_position Legend position. Default is top. One of c("none", "left", "right", "bottom",
#' "top").
#'
#' @export
#'
#' @examples
#' track_reads_dada2(demo_dada2_result$reads_track, single_end = FALSE)

track_reads_dada2 <- function(
  reads_track,
  single_end = FALSE,
  relative_abundance = FALSE,
  legend_position = "top"
  ) {
  if (relative_abundance) {
    reads_track <- reads_track %>% apply(1, function(x) x/x[1]) %>% t()
  }
  reads_track <- rownames_to_column(as.data.frame(reads_track), var = "SampleID")
  reads_track$SampleID <- factor(reads_track$SampleID)
  reads_track <- gather(reads_track, 2:ncol(reads_track), key = "stages", value = "reads")
  if (single_end) {
    reads_track$stages <- factor(
      reads_track$stages, levels = c("input", "filtered", "dereplicated", "nonchim")
      )
  } else {
    reads_track$stages <- factor(
      reads_track$stages, levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
      )
  }
    ggplot(reads_track, aes(x = stages, y = reads, color = SampleID)) +
      scale_color_manual(values = distinctive_colors) +
      geom_point() +
      geom_line(aes(group = SampleID)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 8),
            legend.position = legend_position)
}

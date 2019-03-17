#' dada2_reads_track
#'
#' dada2_reads_track can track the reads count in a dada2 workflow result which created by Xbiome 16S
#' pipeline. Xbiome 16S pipeline dada2 workflow will generate a list that contain sequence table,
#' taxonomy table and reads track data frame. Input the reads track data frame and read type, this
#' function can draw a line plot of reads track of every sample. X-axis will be every stage in dada2
#' workflow, Y-axis will be the reads counts.
#'
#' @param track The reads track data frame from Xbiome 16S pipeline dada2 workflow result.
#' @param single_end Default is FALSE. If single_end == TRUE, means the sequence files are single end,
#'                   the x-axis will contain 'input', 'filtered', 'dereplicated', 'nonchim'. If
#'                   single_end == FALSE, means the sequence files are paired end, the x-axis will
#'                   contain 'input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim'.
#' @export
#' @examples
#' dada2_reads_track(Shaoyifu_dada2_result$reads_track, single_end = FALSE)

dada2_reads_track <- function(track, single_end = FALSE) {
  track <- rownames_to_column(as.data.frame(track), var = "subject_id")
  track$subject_id <- factor(track$subject_id)
  track <- gather(track, 2:ncol(track), key = "stages", value = "reads")
  if (single_end) {
    track$stages <- factor(track$stages, levels = c("input", "filtered",
                                                    "dereplicated", "nonchim"))
  } else {
    track$stages <- factor(track$stages, levels = c("input", "filtered",
                                                    "denoisedF", "denoisedR",
                                                    "merged", "nonchim"))
  }
  ggplot(track, aes(x = stages, y = reads, color = subject_id)) +
    scale_color_manual(values = distinctive_colors) +
    geom_point() +
    geom_line(aes(group = subject_id)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          strip.text.x = element_text(size = 14), 
          # move legend position to the top
          legend.position = "top"
          )
}

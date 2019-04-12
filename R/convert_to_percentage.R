#' convert_to_percentage
#'
#' convert_to_percentage can convert a data frame to percentage.
#'
#' @param df a input data frame.
#' @param row_sum Default is TRUE. If row_sum == TRUE, then will take every value of a row and divide
#'                by the summary of this row, and apply this to every row. If row_sum == FALSE, then
#'                will take every value of a column and divide by the summary of this column, and apply
#'                this to every column.
#' @export
#' @examples
#' convert_to_percentage(demo_dada2_result$seq_tab, row_sum = TRUE) %>% .[,1:5]

convert_to_percentage <- function (df, row_sum = TRUE) {
  # df must be data frame only contain numeric
  if (row_sum) {
    # 将df中每一行求和，然后用每一行中的每一个数字除以行和
    sweep(df, 1, rowSums(df), '/')
  } else {
    # 将df中每一列求和，然后用每一列中的每一个数字除以列和
    sweep(df, 2, colSums(df), '/')
  }
}

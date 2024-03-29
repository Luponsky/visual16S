% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_to_percentage.R
\name{convert_to_percentage}
\alias{convert_to_percentage}
\title{convert_to_percentage}
\usage{
convert_to_percentage(df, row_sum = TRUE)
}
\arguments{
\item{df}{A data frame.}

\item{row_sum}{Default is TRUE. If row_sum == TRUE, then will take every value of a row and divide by
the summary of this row, and apply this to every row. If row_sum == FALSE, then will take every value
of a column and divide by the summary of this column, and apply this to every column.}
}
\description{
convert_to_percentage can convert a data frame to percentage.
}
\examples{
convert_to_percentage(demo_dada2_result$seq_tab, row_sum = TRUE) \%>\% .[,1:5]
}

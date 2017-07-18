#' Process Genome Studio export
#'
#' Processes a Genome Studio data export table to only import useful information for clustering (theta and R values) and ignores the rest
#'
#' @param exportpath the path to the GS data export
#'
#' @return A list of two data.frames. The first is a table of theta values and the second is a table of intensity values


processGSexport <- function(exportpath = NULL) {
  GSheader <-
    as.character(fread(
      exportpath,
      header = F,
      nrows = 1,
      stringsAsFactors = F
    ))
  theta_ind <- grep('\\.Theta', GSheader)
  r_ind <- grep('\\.R', GSheader)
  theta_table <-
    fread(
      exportpath,
      header = T,
      nrows = -1,
      stringsAsFactors = F,
      select = c(2, theta_ind)
    )
  r_table <-
    fread(
      exportpath,
      header = T,
      nrows = -1,
      stringsAsFactors = F,
      select = c(2, r_ind)
    )
  return(list(data.frame(theta_table), data.frame(r_table)))
}

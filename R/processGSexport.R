#' Process Genome Studio export
#'
#' Processes a Genome Studio data export table to only import useful information for clustering (theta and R values) and ignores the rest
#'
#' @param exportpath the path to the GS data export
#'
#' @param transform How to transform theta to control variance. (no transformation = 'none', logit transformation = 'logit', arcsine transformation (default) = 'arcsin')
#'
#' @return A list of two data.frames. The first is a table of theta values and the second is a table of intensity values


processGSexport <- function(exportpath = NULL, transform = 'arcsin') {
  # Read in the header
  GSheader <-
    as.character(fread(
      exportpath,
      header = F,
      nrows = 1,
      stringsAsFactors = F
    ))

  # Read in theta values
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

  # Read in R values
  r_table <-
    fread(
      exportpath,
      header = T,
      nrows = -1,
      stringsAsFactors = F,
      select = c(2, r_ind)
    )

  # Transformation functions
  asinTransform <- function(p) { asin(sqrt(p)) }
  logitTransform <- function(p) { log(p/(1-p)) }

  # Perform theta transformation if specified
  if (transform == 'asin') {theta_table[,2:ncol(theta_table)] <- asinTransform(theta_table[,2:ncol(theta_table)])}
  if (transform == 'logit') {theta_table[,2:ncol(theta_table)] <- logitTransform(theta_table[,2:ncol(theta_table)])}

  # Return both theta and R tables in a list
  return(list(data.frame(theta_table), data.frame(r_table)))
}

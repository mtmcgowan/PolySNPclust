#' extract.clust.index.R
#'
#' extracts the MixmodResults for a given marker index
#'
#' @param index The index of the marker to be extracted
#' @param mixmodout Rmixmod results
#'
#' @return A list of Rmixmod cluster statistics for a particular marker
#'
#' @examples
#'
#'
# This function extracts the MixmodResults for a given marker index
extract.clust.index <- function(index, mixmodout) {
  clust <- mixmodout[[index]][[1]]
  return(clust)
}

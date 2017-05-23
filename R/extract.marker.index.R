#' extract.marker.index
#'
#' extracts the theta and r coordinates for a given marker index
#'
#' @param index The index of the marker to be extracted
#' @param GSdata Processed GenomeStudio data
#'
#' @return A list containing only marker name, Theta and R values for the specified marker
#'
#' @examples
#'
#'
# This function
extract.marker.index <- function(index, GSdata) {
  theta <- unlist(GSdata[[1]][index, -1])
  r <- unlist(GSdata[[2]][index, -1])
  marker_name <- unlist(GSdata[[1]][index, 1])
  marker <- list(marker_name, theta, r)
  return(marker)
}

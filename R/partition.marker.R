#' partition.marker
#'
#' Checks for any missing values removed before clustering, inserts NA values in the mixmod partitioning, and consolidates marker data together for plotting
#'
#' @param index The index of the marker to be extracted
#' @param GSdata Processed GenomeStudio data
#' @param mixmodout Rmixmod results
#'
#' @return A list of data for a specific marker: name, theta, r, partition
#'
#' @examples
#'
#'

partition.marker <- function(index, GSdata, mixmodout) {
  # Extracting the raw results
  marker <- extract.marker.index(index, GSdata)
  clust <- extract.clust.index(index, mixmodout)

  # Extract the cluster assignments
  partition <- clust['partition']

  # Test to see whether any NA values were removed before clustering and re-insert them to make the cluster vector line up with the test_snp matrix
  if (length(mixmodout[[index]]) > 2) {
    nrem <- length(mixmodout[[index]])
    removedsamp <- unlist(mixmodout[[index]][3:nrem])
    for (n in removedsamp) {
      partition <- append(partition, NA, after = n-1)
    }
  }

  # Add the fixed partition to tthe x/y matrix
  marker$cluster <- partition

  names(marker) <- c('Name', 'theta', 'r', 'partition')
  # Return the final matrix
  return(marker)
}

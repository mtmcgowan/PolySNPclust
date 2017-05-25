#' extract.clust.stats
#'
#' Creates a data.frame of cluster statistics for a set of marker indices
#'
#' @param indices A vector of marker indices to extract
#' @param GSdata Processed GenomeStudio data
#' @param mixmodout Clustering output from step 2
#'
#' @return A folder containing a partially populated .csv file and a nested folder of marker plots corresponding to the .csv file
#'
#' @examples
#'

extract.clust.stats <- function(indices, GSdata, mixmodout) {
  # Determine the proper size of the final frame
  clustnum <- unlist(lapply(mixmodout, function(x) {
    clust_stats <- unlist(x[[1]][1])
    return(clust_stats)
  }))

  clustcount <- sum(clustnum)

  # Create the empty matrix
  clust_stats <- data.frame(matrix(NA, nrow = clustcount, ncol = 9))
  names(clust_stats) <- c('marker', 'cluster', 'prop', 'xcoord', 'ycoord', 'var1', 'var2', 'xvar', 'yvar')

  pblapply(indices, function(x) {
    marker <- partition.marker(x, GSdata, mixmodout)
    nc <- mixmodout[[x]][[1]][1]
    clust_temp <- data.frame(matrix(NA, nrow = nc, ncol = 7))
    names(clust_temp) <- c('marker', 'cluster', 'prop', 'xcoord', 'ycoord', 'var1', 'var2')
    model <- mixmodout[[x]][[1]]
    clust_temp$marker <- mixmodout[[x]][[2]]

    for (n in 1:nc) {
      clust_temp$cluster[n] <- n
      # Store cluster number and x/y coordinates
      clust_temp$prop[n] <- model['parameters'][1][n]
      clust_temp$xcoord[n] <- model['parameters'][2][n,1]
      clust_temp$ycoord[n] <- model['parameters'][2][n,2]

      # extract eigen values from the rotational covariance matrix
      eigen_stats <- eigen(model['parameters'][3][[n]])[[1]]
      clust_temp$var1[n] <- eigen_stats[1]
      clust_temp$var2[n] <- eigen_stats[2]

      clust_theta <- marker$theta[marker$cluster == n]
      clust_r <- marker$r[marker$cluster == n]

    }
    return(clust_temp)
  })
}

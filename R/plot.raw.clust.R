#' @title extract.clust.index
#'
#' @description partitions a marker and plots points colored by clusters and surrounded by variance representative elipsoids and cluster centers
#'
#' @param index The index of the marker to be extracted
#' @param GSdata Processed GenomeStudio data, 2 item list of 'theta' and 'r' data frames
#' @param mixmodout Rmixmod results
#'
#' @return A ggplot scatterplot. Colors indicate cluster membership. Red lines indicate cluster center and relative rotation and variance.
#'
#' @examples
#'
#'

plot.raw.clust <- function(index, GSdata, mixmodout) {
  # Extract marker name
  marker_name <- GSdata[[1]][[index, 1]]

  # Extract cluster information
  marker <- partition.marker(index, GSdata, mixmodout)

  # Extract best model
  clust <- extract.clust.index(index, mixmodout)

  # The ggplot scatter plot base
  marker_frame <- data.frame(marker$theta, marker$r, marker$partition)
  names(marker_frame) <- c('theta', 'r', 'cluster')
  bplot <- ggplot(marker_frame, aes(x = theta, y = r, color = factor(cluster)))

  # Calculating elipse data
  ellipse_data <- matrix(nrow = 0, ncol = 3)
  center_data <- matrix(nrow = 0, ncol = 3)
  for (i in 1:clust['nbCluster']) {
    ctr <- clust['parameters'][2][i,]
    A <- clust['parameters'][3][[i]]
    ell_points <- data.frame(car::ellipse(ctr, shape=A, radius=0.98, col="red", lty=2, draw = F))
    ell_points$cluster <- i
    center <- data.frame(t(ctr))
    center$cluster <- i
    ellipse_data <- rbind(ellipse_data, ell_points)
    center_data <- rbind(center_data, center)
  }

  # Consolidating ellipse info into a better format
  ellipse_data <- data.frame(ellipse_data)
  center_data <- data.frame(center_data)
  names(ellipse_data) <- c('x', 'y', 'clust')
  names(center_data) <- c('x', 'y', 'clust')

  # Setting up R colors
  myColors <- brewer.pal(8,"Set2")
  names(myColors) <- levels(1:8)
  colScale <- scale_colour_manual(name = "grp",values = myColors)

  # Plotting the marker
  plot <- bplot + geom_point(size = 0.75) +
    geom_path(data = ellipse_data, aes(x = x, y = y, group = clust), color = 'red', size = 0.5, inherit.aes = F) +
    geom_point(data = center_data, aes(x = x, y = y), shape = 3, color = 'red', size = 2, stroke = 1, inherit.aes = F) +
    colScale +
    ggtitle(paste(marker_name, '\t Index=', index))
  return(plot)
}

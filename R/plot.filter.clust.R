#' @title plot.filter.clust
#'
#' @description Partitions a marker, classifies outlier clusters with a gmb model and plots points colored by clusters and surrounded by variance representative elipsoids and cluster centers
#'
#' @param index The index of the marker to be extracted
#' @param GSdata Processed GenomeStudio data, 2 item list of 'theta' and 'r' data frames
#' @param mixmodout Rmixmod results
#' @param gbm_model A gbm model output from the train.noise.filter command
#'
#' @return A ggplot scatterplot. Colors indicate cluster membership. Red lines indicate cluster center and relative rotation and variance.
#'
#' @examples
#'
#'

plot.filter.clust <- function(index, GSdata, mixmodout, gbm_model) {
  # Extract marker name
  marker_name <- GSdata[[1]][[index, 1]]

  # Extract cluster information
  marker <- partition.marker(index, GSdata, mixmodout)


  # Extract best model
  clust <- extract.clust.index(index, mixmodout)

  # Combine info into a single frame
  marker_frame <- data.frame(marker$theta, marker$r, marker$partition)
  names(marker_frame) <- c('theta', 'r', 'cluster')

  # Extract cluster statistics
  clust_stats <- rbindlist(extract.clust.stats(index, GSdata, mixmodout))

  # Classify outlier clusters
  predictor_names <- c('prop', 'xcoord', 'ycoord', 'var1', 'var2')
  clust_predict <- predict(object=gbm_model, clust_stats[,..predictor_names], type = 'raw')
  clust_stats$class <- clust_predict
  clust_stats$cat_bin <- as.logical(as.numeric(clust_predict)-1)

  # A vector of clusters to keep
  clust_keep <- which(clust_stats$cat_bin)


  # Modify the partition to remove outlier clusters
  filter_partition <- marker$partition
  filter_partition[which(!(filter_partition %in% clust_keep))] <- NA
  removed_points <- marker_frame[is.na(filter_partition) & !is.na(marker_frame$theta),]

  # The ggplot scatter plot base
  bplot <- ggplot(marker_frame, aes(x = theta, y = r, color = factor(cluster)))

  # Calculate max and min values observed
  theta_min <- min(marker_frame$theta, na.rm = T)
  theta_max <- max(marker_frame$theta, na.rm = T)
  r_min <- min(marker_frame$r, na.rm = T)
  r_max <- max(marker_frame$r, na.rm = T)




  # Calculating elipse data
  ellipse_data <- matrix(nrow = 0, ncol = 3)
  center_data <- matrix(nrow = 0, ncol = 3)
  for (i in clust_keep) {
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
  names(myColors) <- levels(c(NA,1:8))
  colScale <- scale_colour_manual(name = "grp",values = myColors)

  # Plotting the marker
  plot <- bplot + geom_point(size = 0.75) +
    geom_path(data = ellipse_data, aes(x = x, y = y, group = clust), color = 'red', size = 0.2, inherit.aes = F) +
    geom_point(data = center_data, aes(x = x, y = y), shape = c(center_data[,3] + 48), color = 'red', size = 3, stroke = 1, inherit.aes = F) +
    colScale +
    ggtitle(paste(marker_name, '\t Index=', index)) +
    geom_point(data = removed_points, aes(x = theta, y = r), color = 'darkgrey', size = 0.75)
  return(plot)
}

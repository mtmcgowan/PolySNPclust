#' plot.markers
#'
#' Plots clustering results for all markers in an index
#'
#' @param index indices of the markers to be plotted
#' @param GSdata Processed GenomeStudio data
#' @param mixmodout Clustering output from step 2
#' @param filter Should all clusters be plotted (pre-filtering), or should the noise filter be applied first. (default = F)
#' @param gbm_model The gbm model to be used for filtering (default = NULL)
#'
#' @return A folder containing plots for all markers in .png form

plot.markers <- function(index, GSdata, mixmodout, filter = F, gbm_model = NULL) {
  cat('Creating plot folder')
  # Make a directory for storing the excel file
  dir.create('marker.plots')

  cat(paste('\n', 'Plotting markers', sep = ''))
  # Plotting the markers
  for (n in index) {
    marker_name <- mixmodout[[n]][[2]]
    file_name <- paste(which(index == n), '_', marker_name, '.png', sep = "")

    if(filter == F) {clustplot <- plot.raw.clust(n, GSdata, mixmodout)}

    if(filter == T) {clustplot <- plot.filter.clust(n, GSdata, mixmodout, gbm_model)}
    ggsave(paste('marker.plots/', file_name, '.png', sep = ""), plot = clustplot, device = 'png')
  }

  cat(paste('\n', 'Finished plotting marker clustering', sep = ''))
}

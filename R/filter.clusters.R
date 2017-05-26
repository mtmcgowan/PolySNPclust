#' Classify unwanted clusters, remove bad clusters, and output partition table
#'
#' Uses a gbm model to predict whether to keep or remove clusters
#'
#' @param filepath the path to the completed .csv file
#'
#' @return A gbm model object
#' @examples
#'

filter.clusters <- function(GSdata, mixmodout, model_gbm) {
  all_clust_stats <- extract.clust.stats(indices = 1:length(mixmodout), GSdata, mixmodout)

  predictor_names <- c('prop', 'xcoord', 'ycoord', 'var1', 'var2')

  clustnum <- unlist(lapply(mixmodout, function(x) {
    clust_stats <- unlist(x[[1]][1])
    return(clust_stats)
  }))

  # Row bind all list items into one large data.table
  clust_stats <- rbindlist(all_clust_stats)

  # Classify clusters
  clust_predict <- predict(object=model_gbm, clust_stats[,..predictor_names], type = 'raw')
  clust_stats$class <- clust_predict
  clust_stats$cat_bin <- as.logical(as.numeric(clust_predict)-1)

  clustnum_filt <- vector(mode = 'integer', length = length(mixmodout))
  for (i in 1:length(mixmodout)) {
    cat('\r', 'Processing marker ', i, ' out of ', length(mixmodout), sep = '')
    marker_name <- mixmodout[[i]][[2]]
    real_clust <- which(clust_stats$marker == marker_name & clust_stats$cat_bin)
    clustnum_filt[i] <- length(real_clust)
  }
  clusthist <- qplot(clustnum_filt, bins = length(unique(clustnum_filt)))
  ggsave(paste('clustnum_hist', '.png', sep = ""), plot = clusthist, device = 'png')
}

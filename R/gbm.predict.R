#' Train noise filter
#'
#' Reads in completed training data and trains a gradient boosted model outputting prediction accuracy statistics
#'
#' @param model_gbm A gradient boosted model object created using manually scored training data
#' @param clust_stats A cluster statistic object from the extract.clust.stats function
#'
#' @return A gbm model object
#' @examples
#' 

gbm.predict <- function(model_gbm, clust_stats) {
  predictor_names <- c('prop', 'xcoord', 'ycoord', 'var1', 'var2')
  preds <- predict(object=model_gbm, clust_stats[,..predictor_names], type = 'raw')
  clust_stats$filter <- preds
  return(clust_stats)
}
#' make.training.data
#'
#' Creates a set of training data for manual scoring using excel from a processed Genome Studio export and an Rmixmod cluster output
#'
#' @param GSdata Processed GenomeStudio data
#' @param mixmodout Clustering output from step 2
#' @param marknum The number of markers to sample for each cluster category (total = 4 * marknum)
#' @param seed Fix the random seed for reproducibility
#' @return A folder containing a partially populated .csv file and a nested folder of marker plots corresponding to the .csv file
#'
#' @examples
#'

make.training.data <- function(GSdata, mixmodout, marknum = 100, seed = NULL) {
  cat('Creating training data folders')
  # Make a directory for storing the excel file
  dir.create('training.data')
  dir.create('training.data/plots')

  setwd('training.data')

  cat('\\n Sampling representative training data markers')
  # Calculate how many clusters were identified for each marker
  clustnum <- unlist(lapply(mixmodout, function(x) {
    clust_stats <- unlist(x[[1]][1])
    return(clust_stats)
  }))

  # Indexing marker category
  clust_2 <- which(clustnum == 2)
  clust_3 <- which(clustnum == 3)
  clust_4 <- which(clustnum == 4)
  clust_5 <- which(clustnum > 5)

  # Randomly sampling within each category (fix seed if specified)
  if(!is.null(seed)){
    set.seed(seed)
  }
  clust_2_samp <- sample(clust_2, marknum)

  if(!is.null(seed)){
    set.seed(seed)
  }
  clust_3_samp <- sample(clust_3, marknum)

  if(!is.null(seed)){
    set.seed(seed)
  }
  clust_4_samp <- sample(clust_4, marknum)

  if(!is.null(seed)){
    set.seed(seed)
  }
  clust_5_samp <- sample(clust_5, marknum)

  # Combine indices together
  training_index <- c(clust_2_samp, clust_3_samp, clust_4_samp, clust_5_samp)

  # Extract cluster information
  training_data_list <- extract.clust.stats(training_index, GSdata, mixmodout)
  training_data <- data.frame(matrix(nrow = 0, ncol = 7))
  for (i in training_data_list) {
    training_data <- rbind(training_data, i)
  }

  training_data$class <- NA

  cat('\\n Writing .csv table of marker statistics')
  write.table(training_data, file = 'training.csv', row.names = F, quote = F, col.names = T, sep = ",")

  cat('\\n Plotting markers')
  # Plotting the markers
  for (n in training_index) {
    marker_name <- mixmodout[[n]][[2]]
    file_name <- paste(which(training_index == n), '_', marker_name, '.png', sep = "")
    clustplot <- plot.raw.clust(n, GSdata, mixmodout)
    ggsave(paste('plots/', file_name, '.png', sep = ""), plot = clustplot, device = 'png')
  }

  # return to the parent directory
  setwd('..')

  cat('\\n Completed creation of training data form and plots')
}

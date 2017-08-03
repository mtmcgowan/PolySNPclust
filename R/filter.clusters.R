#' Classify unwanted clusters, remove bad clusters, and output partition table
#'
#' Uses a gbm model to predict whether to keep or remove clusters
#'
#' @param filepath the path to the completed .csv file
#'
#' @return A gbm model object
#' @examples
#'

filter.clusters <- function(GSdata, mixmodout, gbm_model) {
  cat('\n', 'Extracting cluster statistics for all markers', '\n')
  all_clust_stats <- extract.clust.stats(indices = 1:length(mixmodout), GSdata, mixmodout)

  predictor_names <- c('prop', 'xcoord', 'ycoord', 'var1', 'var2')

  clustnum <- unlist(lapply(mixmodout, function(x) {
    clust_stats <- unlist(x[[1]][1])
    return(clust_stats)
  }))

  # Row bind all list items into one large data.table
  clust_stats <- rbindlist(all_clust_stats)

  # Classify clusters
  cat('\n', 'Classifying clusters with gbm model', '\n')
  clust_predict <- predict(object=gbm_model, clust_stats[,..predictor_names], type = 'raw')
  clust_stats$class <- clust_predict
  clust_stats$cat_bin <- as.logical(as.numeric(clust_predict)-1)

  # Count the number of clusters that pass the filter for each marker
  clustnum_filt <- unlist(pblapply(1:length(mixmodout), function(x) {
    marker_name <- mixmodout[[x]][[2]]
    real_clust <- which(clust_stats$marker == marker_name & clust_stats$cat_bin)
    return(length(real_clust))
  }))

  # Generate a histogram plot showing the distribution for post-filter cluster count
  cat('\n', 'Plotting and saving the cluster count distribution', '\n')
  clusthist <- qplot(clustnum_filt, bins = length(unique(clustnum_filt)))
  ggsave(paste('clustnum_hist', '.png', sep = ""), plot = clusthist, device = 'png')

  # Bin markers into categories based on cluster count
  cat('\n', 'Seperating markers into cluster count categories', '\n')
  one_clust <- which(clustnum_filt == 1)
  two_clust <- which(clustnum_filt == 2)
  three_clust <- which(clustnum_filt == 3)
  fourplus_clust <- which(clustnum_filt > 3)

  three_clust_names <- as.character(GSdata[[1]]$Name[three_clust])
  three_clust_stats <- clust_stats[clust_stats$marker %in% three_clust_names & clust_stats$cat_bin == T,]

  # Calculate the minimum distance and distance ratio values for three-cluster markers
  three_clust_ratio <- vector(mode = 'integer', length = length(three_clust_names))
  three_clust_mindist <- vector(mode = 'integer', length = length(three_clust_names))

  for (i in 1:length(three_clust_names)) {
    xcoord <- three_clust_stats$xcoord[three_clust_stats$marker == three_clust_names[i]]
    distances <- sort(dist(xcoord))
    ratio <- distances[1] / distances[2]
    three_clust_mindist[i] <- min(distances)
    three_clust_ratio[i] <- ratio
  }

  # Defining a short function for setting cutpoints
  define_cp <- function(obj, xname) {
    dist_hist <- hist(obj, breaks = round(length(obj) * 0.05), main = paste("Histogram of" , xname) )
    coords <- locator(type="l", n = 1)
    return(coords$x)
  }

  # Have the user set the mindist cutpoint for three-cluster markers
  cat('\n', 'Click where to set the minimum distance cutpoint for three-cluster markers', '\n')
  mindist_cp <- define_cp(three_clust_mindist, xname = 'Minimum cluster distance')

  # Have the user set the ratio cutpoint for three-cluster markers
  cat('\n', 'Click where to set the minimum ratio cutpoint for three-cluster markers', '\n')
  ratio_cp <- define_cp(three_clust_ratio, xname = 'Cluster distance ratio')

  # Calculate the minimum distance and distance ratio values for three-cluster markers
  two_clust_names <- as.character(GSdata[[1]]$Name[two_clust])
  two_clust_stats <- clust_stats[clust_stats$marker %in% two_clust_names & clust_stats$cat_bin == T,]

  two_clust_mindist <- vector(mode = 'integer', length = length(two_clust_names))

  for (i in 1:length(two_clust_names)) {
    xcoord <- two_clust_stats$xcoord[two_clust_stats$marker == two_clust_names[i]]
    distances <- sort(dist(xcoord))
    two_clust_mindist[i] <- min(distances)
  }

  # Have the user set the mindist cutpoint for two-cluster markers
  cat('\n', 'Click where to set the minimum ratio cutpoint', '\n')
  mindist_cp2 <- define_cp(two_clust_mindist, xname = 'Minimum cluster distance')

  cat('\n', 'Applying filters', '\n')

  # Use defined cutpoints to filter the markers
  three_clust_filter_list <- three_clust_names[which(three_clust_mindist >= mindist_cp & three_clust_ratio >= ratio_cp)]
  two_clust_filter_list <- two_clust_names[which(two_clust_mindist >= mindist_cp2)]

  # Combine both vectors into a single list
  markerlist <- c(two_clust_filter_list, three_clust_filter_list)

  # Calling the extract.genotypes function to create the genotype matrix
  genotype_matrix <- extract.genotype(markerlist, GSdata, mixmodout, clust_stats)

  # Adding sample names to the genotype matrix
  genotype_matrix$taxa <- names(GSdata[[1]])[-1]
  genotype_matrix <- genotype_matrix[, c(ncol(genotype_matrix), 1:ncol(genotype_matrix)-1)]

  return(list(genotype_matrix, clust_stats))
}


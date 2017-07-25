#' Cluster markers in Parallel with Rmixmod
#'
#' Sets up a parallel backend and clusters each marker using mixture of gaussian models provided by the 'Rmixmod' package.
#'
#' @param GSdata Processed GenomeStudio data, 2 item list of 'theta' and 'r' data frames
#'
#' @return An Rmixmod object containing the best cluster model results
#'
#' @example mixmodout <- cluster.markers(GSdata)


cluster.markers <- function(GSdata) {

  # Convert GS list to individual tables (to make clustering functions simpler)
  theta_frame <- GSdata[[1]]
  r_frame <- GSdata[[2]]

  rownames(theta_frame) <- unlist(GSdata[[1]][,1])
  rownames(r_frame) <- unlist(GSdata[[2]][,1])

  # Extracting the list of markers to iterate over
  marker_list <- theta_frame$Name

  # Initializing sub-functions
  # A function for extracting marker information
  # Parameters: x (marker name), theta_frame (matrix contining theta values), r_frame (matrix containing r values).
  # Output: a matrix of x and y vectors (continins only theta and r values for that particular marker)
  snp_extract <- function(x, theta_frame, r_frame) {
    w <- which(theta_frame$Name == x)
    theta <- unlist(theta_frame[w,2:ncol(theta_frame)])
    r <- unlist(r_frame[w,2:ncol(theta_frame)])
    test_snp <- data.frame(theta, r)
    names(test_snp) <- c('x', 'y')
    return(test_snp)
  }

  # A function that will take a test_snp consisting of theta and R values and cluster using Rmixmod
  # 3 parameters: test_snp (output from snp_extract), model_list (a list of gaussian models to run), clustnum = (vector of # of clusters to test)
  # Output: a list of MixmodCluster objects for a particular marker (includes sub-optimal models)
  snp_mixclust <- function(test_snp, model_list = c("Gaussian_pk_L_Ck", "Gaussian_pk_L_Bk"), clustnum = 1:8) {

    clust_strategy <- mixmodStrategy(algo = c('EM', 'CEM'), nbTry = 20,
                                     initMethod = "CEM",
                                     nbTryInInit = 1000, nbIterationInInit = 5,
                                     nbIterationInAlgo = 500, epsilonInInit = 0.001,
                                     epsilonInAlgo = 0.001, seed = NULL, parameter=NA,
                                     labels=NA)
    clust_results <- mixmodCluster(test_snp, nbCluster = clustnum, models = mixmodGaussianModel(listModels = model_list), criterion = 'ICL', strategy = clust_strategy)
    return(clust_results)
  }

  # A function that combines the above two functions and extracts only the 'best
  # Parameters: x (the marker name), theta_frame (the data frame containing theta values), R_frame (the data frame containing R values)
  # Output: A MixmodCluster list of data for the best cluster model
  marker_clust <- function(x) {
    place <- which(marker_list == x)
    total <- length(marker_list)
    print(paste('Processing ', place, ' out of ', total))

    y <- snp_extract(x, theta_frame, r_frame)

    # Check for outliers using a box-plot strategy (anything outlide 1.5*IQR)
    bplot <- boxplot.stats(y$y)
    outliers <- which(y$y %in% bplot$out)
    outliers_names <- row.names(y)[outliers]

    # Identify NA values
    NA_samples <- row.names(y[is.na(y$x) | is.na(y$y),])

    # Create a list of samples to remove
    bad_samples <- c(outliers_names, NA_samples)
    bad_samp_num <- length(bad_samples)

    # Also create a vector of which indices are bad values
    samp_rm <- which(row.names(y) %in% bad_samples)

    # Remove bad samples
    if (bad_samp_num > 0) {test_snp <- y[-which(row.names(y) %in% bad_samples),]} else {test_snp <- y}

    clust_results <- snp_mixclust(test_snp)
    model <- clust_results['bestResult']
    marker_name <- x
    b <- list(model, marker_name, samp_rm)
    return(b)
  }

  # Setting a 'cl' variable to the number of logic units available
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, outfile = '')

  # Exporting the required functions and data.frames to each logic unit before running the parallel function
  clusterExport(cl, list('theta_frame', 'r_frame', 'marker_list', "snp_extract", "marker_clust", "snp_mixclust", "mixmodCluster", 'mixmodGaussianModel', 'mixmodStrategy'), envir=environment())
  ptm <- proc.time()
  parallel_out <- parLapply(
    cl, marker_list, marker_clust
  )
  parallel_21_time <- proc.time() - ptm

  stopCluster(cl)

  return(parallel_out)

}

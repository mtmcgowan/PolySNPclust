#' Generates a c(0,1,2) genotype matrix from a list of marker names
#'
#' Uses a gbm model to predict whether to keep or remove clusters
#'
#' @param markerlist A character vector of marker names
#' @param GSdata A list of R and Theta data.frames
#' @param mixmodout Rmixmod clustering results
#' @param clust_stats A data.frame containing statistics for all clusters
#'
#' @return A numerical matrix with rows = markers and columns = samples


extract.genotype <- function(markerlist, GSdata, mixmodout, clust_stats) {
  # Create an empty matrix to store genotype vectors in
  genotype_matrix <- data.frame(matrix(nrow = ncol(GSdata[[1]]) -1, ncol = (length(markerlist))))

  # Iterate through each marker
  for (i in 1:length(markerlist)) {
    # Print a status update
    cat('\r', 'Processing marker ', i, ' out of ', length(markerlist), sep = '')
    # Store marker name
    marker_name <- markerlist[i]
    names(genotype_matrix)[i] <- marker_name
    index <- which(GSdata[[1]]$Name == marker_name)

    # Extract unfiltered partition
    marker <- partition.marker(index, GSdata, mixmodout)

    # Extract best model
    clust <- extract.clust.index(index, mixmodout)

    # Determining which clusters were filtered out and replace bad clusters with NA
    clust_filt <- clust_stats[clust_stats$marker == marker_name,]
    clust_keep <- which(clust_filt$cat_bin)
    marker$partition[!marker$partition %in% clust_keep] <- NA

    clustnum <- sum(!is.na(unique(marker$partition)))

    if (clustnum == 2) {
      good_clust <- clust_filt[clust_filt$cat_bin,]
      clust_order <- good_clust[order(good_clust$xcoord),]
      marker$genotype[marker$partition == clust_order$cluster[1]] <- 0
      marker$genotype[marker$partition == clust_order$cluster[2]] <- 2
      genotype_matrix[,i] <- marker$genotype
    }
    else if (clustnum == 3) {
      good_clust <- clust_filt[clust_filt$cat_bin,]
      clust_order <- good_clust[order(good_clust$xcoord),]
      marker$genotype[marker$partition == clust_order$cluster[1]] <- 0
      marker$genotype[marker$partition == clust_order$cluster[2]] <- 1
      marker$genotype[marker$partition == clust_order$cluster[3]] <- 2
      genotype_matrix[,i] <- marker$genotype
    }
  }
  return(genotype_matrix)
}

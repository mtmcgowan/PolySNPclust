#' scrub.GSdata
#' @title scrub.GSdata
#'
#' @description Checks for problematic samples and markers. Samples with a set missing rate are flagged. Markers with zero-variance (either all theta or R values are the same) are also flagged for removal.
#'
#' @param GSdata Processed GenomeStudio data, 2 item list of 'theta' and 'r' data frames
#'
#' @param missrate The missing rate threshold where markers or samples will be removed (default = 0.2 i.e. 20% cutoff)
#' @return A item list of 'theta' and 'r' data frames with flagged markers/samples removed.
#'
#' @example GSdata_clean <- scrub.GSdata(GSdata)


scrub.GSdata <- function(GSdata, missrate = 0.2) {
  # Calculate 'theta' and 'r' missing rates for SAMPLES
  samp_theta_missrate <- apply(GSdata[[1]][,2:ncol(GSdata[[1]])], 2, function(x) {sum(is.na(x)) / length(x)})
  samp_r_missrate <- apply(GSdata[[2]][,2:ncol(GSdata[[2]])], 2, function(x) {sum(is.na(x)) / length(x)})

  # Determine which SAMPLES should be removed
  theta_samp_rm <- which(samp_r_missrate > missrate)
  r_samp_rm <- which(samp_r_missrate > missrate)
  samp_rm <- unique(c(theta_samp_rm, r_samp_rm))
  samp_rm_names <- names(GSdata[[1]])[samp_rm]

  # Calculate 'theta' and 'r' missing rates for MARKERS
  mark_theta_missrate <- apply(GSdata[[1]][,2:ncol(GSdata[[1]])], 1, function(x) {sum(is.na(x)) / length(x)})
  mark_r_missrate <- apply(GSdata[[2]][,2:ncol(GSdata[[2]])], 1, function(x) {sum(is.na(x)) / length(x)})

  # Calculate variance for 'theta' and 'r' variables for MARKERS
  mark_theta_var <- apply(GSdata[[1]][,2:ncol(GSdata[[1]])], 1, function(x) {var(x, na.rm = T)})
  mark_r_var <- apply(GSdata[[2]][,2:ncol(GSdata[[2]])], 1, function(x) {var(x, na.rm = T)})

  # Determine which MARKERS should be removed
  mark_theta_rm <- which(mark_r_missrate > missrate | mark_theta_var == 0)
  mark_r_rm <- which(mark_r_missrate > missrate | mark_r_var == 0)
  mark_rm <- unique(c(mark_theta_rm, mark_r_rm))
  mark_rm_names <- (GSdata[[1]][mark_rm, 1])

  # Subset GSdata by removing bad samples and markers
  if (length(samp_rm) == 0 & length(mark_rm) == 0) {
    theta_clean <- GSdata[[1]]
    r_clean <- GSdata[[2]]
    GSdata_clean <- list(theta_clean, r_clean)
    cat("No samples or markers removed", sep = "\n")
  }


  else if(length(samp_rm) == 0 & length(mark_rm) > 0) {
    theta_clean <- GSdata[[1]][-c(mark_rm),]
    r_clean <- GSdata[[2]][-c(mark_rm),]
    GSdata_clean <- list(theta_clean, r_clean)
    cat("No samples removed", sep = "\n")
    cat("Removing markers:", c(mark_rm_names), sep = "\n")
  }

  else if(length(samp_rm) > 0 & length(mark_rm) == 0) {
    theta_clean <- GSdata[[1]][,-c(samp_rm+1)]
    r_clean <- GSdata[[2]][,-c(samp_rm+1)]
    GSdata_clean <- list(theta_clean, r_clean)
    cat("Removing samples:", c(samp_rm_names), sep = "\n")
    cat("No markers removed", sep = "\n")
  }

  else if(length(samp_rm) > 0 & length(mark_rm) > 0) {
    theta_clean <- GSdata[[1]][-c(mark_rm),-c(samp_rm+1)]
    r_clean <- GSdata[[2]][-c(mark_rm),-c(samp_rm+1)]
    GSdata_clean <- list(theta_clean, r_clean)
    cat("Removing markers:", c(mark_rm_names), sep = "\n")
    cat("Removing samples:", c(samp_rm_names), sep = "\n")
  }

  return(GSdata_clean)
}






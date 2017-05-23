#' SLURMclustpar
#'
#' Create a single SLURM script for submitting parallel jobs across multiple Kamiak nodes
#'
#' @param GSdata Processed GenomeStudio data ready for clustering
#' @param partition The partition to use on Kamiak
#' @param account The account to use on Kamiak
#' @param nodes The number of nodes to use
#' @param cpus The number of cpus per node to use
#' @param job The job ID to use for file creation
#' @return A folder containing the SLURM script and sub-directories containing tables for clustering.
#'
#' @examples

SLURM.clust.par <- function(GSdata, partition, account, nodes, cpus, job) {
  dir.create(job)
  set.wd(job)

  for (i in 1:nodes) {

  }
}

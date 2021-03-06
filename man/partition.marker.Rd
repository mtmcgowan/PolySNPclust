% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition.marker.R
\name{partition.marker}
\alias{partition.marker}
\title{partition.marker}
\usage{
partition.marker(index, GSdata, mixmodout)
}
\arguments{
\item{index}{The index of the marker to be extracted}

\item{GSdata}{Processed GenomeStudio data}

\item{mixmodout}{Rmixmod results}
}
\value{
A list of data for a specific marker: name, theta, r, partition
}
\description{
Checks for any missing values removed before clustering, inserts NA values in the mixmod partitioning, and consolidates marker data together for plotting
}
\examples{


}

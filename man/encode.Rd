% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadShape.R
\name{encode}
\alias{encode}
\title{encode meter data load shapes according to an existing dictionary}
\usage{
encode(rdata, dic, relerror = 0)
}
\arguments{
\item{rdata}{raw data as a matrix with each row as a day of load data}

\item{dic}{dictionary of cluster centers to use for encodings}

\item{relerror}{whether to include relative error information in result}
}
\value{
encoding: closest shape code, daily sum, L2 err on normalized data, estimated threshold (relative L2 error)
}
\description{
given a matrix of raw meter data (one day per row), encodes each day into the closest fit culster
}


# Copyright 2016 The Board of Trustees of the Leland Stanford Junior University.
# Direct inquiries to Sam Borgeson (sborgeson@stanford.edu)
# or professor Ram Rajagopal (ramr@stanford.edu)

#' Smart Meter load shape data analysis tools for R
#'
#' @section Core functions:
#'
#' VISDOM load shape relies most heavily on the following core functions. See the example vignettes for examples of their usage:
#'
#' \itemize{
#'   \item \code{\link{visdom}}: Depends lightly on the core visdom package functions, mainly the ability for a data source to return data in the appropriate format.
#'
#'
#'   \item \code{\link{create_dictionary}}: Starting with raw load data, create a dictionary of cluster centers and return encodings in one shot.
#'
#'
#'   \item \code{\link{dictionary.distance}}: Compare two load shape dictionaries.
#'
#'
#'   \item \code{\link{draw.top.n.shapes}}: Given a set of encodings and their corresponding dictionary, plot the top load shapes present.
#'
#'
#'   \item \code{\link{encode}}: Given raw meter data and a target dictionary, encode each load day worth of data with its best fit cluster assignment.
#'
#'
#'   \item \code{\link{impute}}: Impute meter data to eliminate blanks (necessary for clustering algorithm to function).
#'
#'
#'   \item \code{\link{quick.norm}}: Efficiently normalize each day of load data for use in a k-mean family algorithm.
#'
#'
#'   \item \code{\link{raw2encoded}}: Starting with raw meter data, cluster into a dictionary of best fit load shapes and perform encoding with that dictionary.
#'
#'
#'   \item \code{\link{reduce.dictionary}}: Shrink a dictionary by merging the most similar entries with each other.
#'
#'
#'   \item \code{\link{shannon.entropy}}: Information entropy calculation on a table of load shape encoding counts.
#'
#'
#'   \item \code{\link{shannon.entropy2}}: Information entropy calculation on load shape encodings.
#'
#' }
#'
#' @docType package
#' @name visdomloadshape
NULL

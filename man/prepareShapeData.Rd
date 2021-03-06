% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadShape.R
\name{prepareShapeData}
\alias{prepareShapeData}
\title{Prepare load shape data for load shape clustering}
\usage{
prepareShapeData(rawData, forceSameDuration = F, subtractMins = F,
  minPower = NULL, coarseTimePeriods = F, return.mins = F)
}
\arguments{
\item{rawData}{data.frame of load shape data, with the first n columns assumed to be metadata, which must include an 'id' column, and the last 24 columns as load shape data.}

\item{forceSameDuration}{applies a heuristic to preserve only data from ids that have the modal number of days of meter data.}

\item{subtractMins}{'de-mins' data by subtracting the daily min observation from all daily observations.}

\item{minPower}{if not null, the function removes load shapes that average less than the provided minPower level.}

\item{coarseTimePeriods}{down samples data by averaging 24 columns of meter data observation to 4 6-hour averages.}

\item{return.mins}{whether or not to return the daily min values subtracted by subtractMins}
}
\description{
Offers several options for cleaning and preparing load shape data, including forcing all data
to be the same duration, subtracting the daily minimum from each day of load, eliminating load shapes that
average less than a minimum power threshold, and down-sampling load data from 24 observations a day to 4.
}


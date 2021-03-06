% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadShape.R
\name{dictionary.distance}
\alias{dictionary.distance}
\title{calculate the distance btw two dictionaries}
\usage{
dictionary.distance(A, B, pa = 1, pb = 1, emd = 1, same = 0,
  tmpdist = NULL)
}
\arguments{
\item{A}{first dictionary: n1 by p matrix}

\item{B}{second dictionary: n2 by p matrix}

\item{pa}{probability distribution vector of dictionary A}

\item{pb}{probability distribution vector of dictionary B}

\item{emd}{whether to estimate the distance btw two load shapes as EMD (earth mover distance) or L1 distance/2}

\item{same}{whether A and B are the same dictionaries or not}

\item{tmpdist}{user can provide the distance matrix (n1 by n2) btw the dictionaries A and B if already calculated}
}
\description{
calculate the distance btw two dictionaries
}


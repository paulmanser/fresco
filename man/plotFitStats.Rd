\name{plotFitStats}
\alias{plotFitStats}
\title{Get fit statistics R^2, AICC, and GCV for loess fits}
\usage{
plotFitStats(object, useControls = TRUE, loessSpan = seq(0.05, 0.95, 0.15),
  sdThreshold = 0.1)
}
\arguments{
  \item{object}{a \code{MethylSet} object}

  \item{useControls}{Should empirical controls be used to
  align and fit loess surfaces?}

  \item{loessSpan}{Supply vector of possible spans for
  fitting loess surface}

  \item{sdThreshold}{Threshold to filter empirical controls
  by standard deviation}
}
\description{
Get fit statistics R^2, AICC, and GCV for loess fits
}


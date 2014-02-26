\name{returnFitStats}
\alias{returnFitStats}
\title{Get fit statistics R^2, AICC, and GCV for loess fits}
\usage{
returnFitStats(object, useControls = TRUE, loessSpan = 0.1,
  sdThreshold = 0.1)
}
\arguments{
  \item{object}{\code{MethylSet} object}

  \item{useControls}{Should empirical controls be used to
  align and fit loess surfaces?}

  \item{loessSpan}{Supply span for fitting loess surface}

  \item{sdThreshold}{Threshold to filter empirical controls
  by standard deviation}
}
\description{
Get fit statistics R^2, AICC, and GCV for loess fits
}


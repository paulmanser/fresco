\name{preprocessFresco}
\alias{preprocessFresco}
\title{Normalize data using flexible local regression using empirical controls}
\usage{
preprocessFresco(object, useControls = TRUE, loessSpan = 0.1,
  fitLoess = TRUE, sdThreshold = 0.1, verbose = TRUE)
}
\arguments{
  \item{object}{\code{MethylSet} object}

  \item{useControls}{Should empirical controls be used to
  align and fit loess surfaces?}

  \item{loessSpan}{Supply span for fitting loess surface}

  \item{fitLoess}{Should loess curve be fitted after
  initial alignment and scaling?}

  \item{sdThreshold}{Threshold to filter empirical controls
  by standard deviation}
}
\description{
Normalize data using flexible local regression using
empirical controls
}


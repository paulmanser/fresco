\name{frescoData}
\alias{frescoData}
\docType{data}
\title{
Data for use with fresco package
}
\description{
Contains an index of empirical control probes used for loess surface fitting.
Contains information on minor allele frenquencies for probes containing SNPS as
well as that cross-hybridize to probes on sex chromosomes obtained from Chen et al.
}
\usage{data(frescoData)}
\format{
  A data frame with 485512 observations on the following 10 variables.
  \describe{
    \item{\code{probeType}}{a factor with levels \code{I} \code{II}}
    \item{\code{targetGC}}{a numeric vector}
    \item{\code{crossHyb}}{a logical vector}
    \item{\code{chromosome}}{a factor with levels \code{} \code{1} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{2} \code{20} \code{21} \code{22} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{X} \code{Y}}
    \item{\code{AF}}{a numeric vector}
    \item{\code{AFR_AF}}{a numeric vector}
    \item{\code{AMR_AF}}{a numeric vector}
    \item{\code{ASN_AF}}{a numeric vector}
    \item{\code{EUR_AF}}{a numeric vector}
    \item{\code{eControls}}{a factor with levels \code{Hemimethylated} \code{Methylated} \code{Unmethylated}}
  }
}
\details{
frescoData is called from within the fresco functions, but a separate copy can be
loaded into the global environment if desired.
}

\examples{
}
\keyword{datasets}

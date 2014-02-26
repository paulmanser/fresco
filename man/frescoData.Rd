\name{frescoData}
\alias{frescoData}
\docType{data}
\title{
Data for use with preprocessFresco
}
\description{
An index of empirical control probes used for loess surface fitting.
Contains information on minor allele frenquencies for probes containing SNPS as
well as probes known to cross-hybridize to sex chromosomes.
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(frescoData)
## maybe str(frescoData) ; plot(frescoData) ...
}
\keyword{datasets}

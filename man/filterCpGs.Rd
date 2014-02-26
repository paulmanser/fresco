\name{filterCpGs}
\alias{filterCpGs}
\title{Filter unreliable probes}
\usage{
filterCpGs(object, removeChromosomes = NULL, filterCrossHyb = TRUE,
  filterNA = TRUE, filterSNP = TRUE, minorAlleleFreq = 0,
  population = "All")
}
\arguments{
  \item{object}{\code{MethylSet} object}

  \item{removeChromosomes}{A character string of chromoses
  to remove}

  \item{filterCrossHyb}{Filter autosomal probes that
  cross-hybridize to sex chromosomes?}

  \item{filterNA}{Filter probes containing at least one
  NA?}

  \item{filterSNP}{Filter probes containing SNPs?}

  \item{minorAlleleFreq}{What is the largest minor allele
  frequency we are willing to tolerate?}

  \item{population}{What population should be used to
  compute minor allele frequency? Default is 'All'}
}
\description{
Filter unreliable probes
}


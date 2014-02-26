#' Filter unreliable probes 
#' 
#' @param object \code{MethylSet} object
#' @param removeChromosomes A character string of chromoses to remove
#' @param filterCrossHyb Filter autosomal probes that cross-hybridize to sex chromosomes?
#' @param filterNA Filter probes containing at least one NA?
#' @param filterSNP Filter probes containing SNPs?
#' @param minorAlleleFreq What is the largest minor allele frequency we are willing to tolerate?
#' @param population What population should be used to compute minor allele frequency?
#' Default is 'All'
#' 
#' @export filterCpGs

filterCpGs <- function(object, removeChromosomes = NULL, filterCrossHyb = TRUE, 
                       filterNA = TRUE, filterSNP = TRUE, 
                       minorAlleleFreq = 0, population = 'All'){
  
  
  if (!is(object, "MethylSet")) stop("'object' needs to be a 'MethylSet'")
  
  populationAF <- c('All', 'African', 'American', 'Asian', 'European')
  
  if (!population %in% populationAF){
    stop("population' must be one of 'All', 'African', 'American', 'Asian', or 'European'")
  }
  
  if (sum(!removeChromosomes %in% c('X', 'Y', 1:22)) > 0){
    stop("'removeChromosomes' needs to be a list of
         chromosomes to remove e.g. c('X', '1')")
  }
  
  data(frescoData)
  removeProbes <- NULL
  probeIDs <- rownames(getMeth(object))
  frescoData <- frescoData[match(probeIDs, rownames(frescoData)), ]
  methSignals <- getMeth(object)
  unmethSignals <- getUnmeth(object)
  
  if (length(removeChromosomes) > 0){
    removeProbes <- c(removeProbes, probeIDs[which(frescoData$chromosome %in% removeChromosomes)])
  }
  
  if (filterCrossHyb){
    removeProbes <- c(removeProbes, probeIDs[which(frescoData$crossHyb)])
  }
  
  if (filterNA){
    NAind <- probeIDs[which(rowSums(is.na(getBeta(object))) > 0)]
    removeProbes <- c(removeProbes, probeIDs[NAind])
  }
  
  if (filterSNP){
    AFtype <- match(population, populationAF) + 4
    SNPind <- which(frescoData[, AFtype] > minorAlleleFreq)
    removeProbes <- c(removeProbes, probeIDs[SNPind])
  }
    
  removeProbes <- unique(removeProbes)
  
  keepCpGs <- setdiff(probeIDs, removeProbes)
  
  out <- object
  methSignals <- methSignals[keepCpGs, ]
  unmethSignals <- methSignals[keepCpGs, ]
  assayDataElement(out, 'Unmeth') <- unmethSignals
  assayDataElement(out, 'Meth') <- methSignals
  
  out  
  
}

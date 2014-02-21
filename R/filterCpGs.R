filterCpGs <- function(object, removeChromosomes = NULL, filterCrossHyb = TRUE, 
                       filterSNP = TRUE, minorAlleleFreq = 0, population = 'All'){
  
  if (!is(object, "MethylSet")) stop("'object' needs to be a 'MethylSet'")
  
  populationAF <- c('All', 'African', 'American', 'Asian', 'European')
  removeProbes <- NULL
  
  probeIDs <- rownames(getMeth(object))
  frescoData <- frescoData[match(probeIDs, rownames(frescoData)), ]
  methSignals <- getMeth(object)
  unmethSignals <- getUnmeth(object)
  
  if (length(removeSexChr) > 0){
    removeProbes <- c(removeProbes, probeIDs[which(frescoData$chromosome %in% removeChr)])
  }
  
  if (filterCrossHyb){
    removeProbes <- c(removeProbes, probeIDs[which(frescoData$crossHyb)])
  }
  
  if (filterSNP){
    AFtype <- match(population, populationAF) + 4
    SNPind <- which(frescoData[,AFtype] > minorAlleleFreq)
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

#' Quality assessment and filtering threshold selection for empirical controls
#' 
#' @param object \code{MethylSet} object
#' @param sdThreshold Standard deviation cut-off for filtering empirical controls
#' 
#' @export empiricalControlQA


empiricalControlQA <- function(object, sdThreshold = .1){

  if (!is(object, "MethylSet")) stop("'object' needs to be a 'MethylSet'")
  
  data(frescoData)
  
  # pull out control probes and get average --------------------------------------
  betaVals <- getBeta(object)
  frescoData <- frescoData[match(rownames(betaVals), rownames(frescoData)), ]
  controlInd <- which(!is.na(frescoData$eControls))
  controlInd <- intersect(controlInd, which(rowSums(is.na(betaVals)) == 0))
  means <- rowMeans(betaVals[controlInd,])  
  
  # plot sorted control probes as heat map -------------------------------------
  par(mfrow = c(1, 3))
  image(betaVals[controlInd[order(means)], ], axes = FALSE,
        main = 'Empirical Control Probe QC',
        xlab = 'CpGs ordered by avg methylation',
        ylab = 'Samples')
  lines(seq(0, 1, length.out = length(means)), 
        means[order(means)], col = 1)
  
  # plot control probe standard deviations --------------------------------------
  controlsSD <- rowSds(betaVals[controlInd[order(means)]])
  
  plot(density(controlsSD), 
       main = 'Empirical Control Probe Standard Deviations',
       xlab='Standard Deviation')
    
  abline(v = sdThreshold)
  legend('right', legend = paste(sum(controlsSD < sdThreshold), 'of', 
                                 length(controlInd), 'controls remaining'), bty = 'n')
  
  image(betaVals[controlInd[order(means)], ][which(controlsSD < sdThreshold), ], 
        axes = FALSE,
        main = 'Filtered Empirical Control Probes',
        xlab = 'CpGs ordered by avg methylation',
        ylab = 'Samples')
  
  lines(seq(0, 1, length.out = length(which(controlsSD < sdThreshold))), 
        means[order(means)][which(controlsSD < sdThreshold)], col = 1)
}

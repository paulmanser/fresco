empiricalControlQA <-
function(object, newControlSet=FALSE, sdCutoff = .1){

  require(minfi)
  if (!is(object, "RGChannelSet")) stop("'object' needs to be a 'RGChannelSet'")
  
  # pull out control probes and get average --------------------------------------
  betaVals <- getBeta(object)
  frescoData <- frescoData[match(rownames(betaVals), rownames(frescoData)), ]
  controlInd <- which(!is.na(frescoData$eControls))
  controlInd <- intersect(controlInd, which(rowSums(is.na(betaVals)) == 0))
  means <- rowMeans(betaVals[controlInd,])  
  
  # plot sorted control probes as heat map -------------------------------------
  image(betaVals[controlInd[order(means)], ], axes = FALSE,
        main = 'Empirical Control Probe QC',
        xlab = 'CpGs ordered by avg methylation',
        ylab = 'Samples')
  lines(seq(0, 1, length.out = length(means)), 
        means[order(means)], col = 1)
  
  # plot control probe standard deviations --------------------------------------
  controlsSD <- apply(betaVals[controlInd[order(means)], ], 1, sd)
  
  plot(density(controlsSD), main = 'Empirical Control Probe Standard Deviations',
       xlab='Standard Deviation')
  
  if (newControlSet){
    
    abline(v=sdCutoff)
    legend('right', legend = paste(sum(controlsSD < sdCutoff), 'of', 
                                   length(controlInd), 'controls remaining'), bty = 'n')
    
    image(betaVals[controlInd[order(means)], ][which(controlsSD < sdCutoff), ], 
          axes=FALSE,
          main='Filtered Empirical Control Probes',
          xlab='CpGs ordered by avg methylation',
          ylab='Samples')
    
    lines(seq(0, 1, length.out = length(which(controlsSD < sdCutoff))), 
          means[order(means)][which(controlsSD < sdCutoff)], col = 1)

    discardControls <- names(controlsSD[controlsSD > sdCutoff])
    frescoData$eControlsFiltered <- frescoData$eControls
    frescoData$eControlsFiltered[rownames(frescoData) %in% discardControls] <- NA
    return(frescoData)

  }
}

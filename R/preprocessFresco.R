preprocessFresco <-
function(object, useControls = TRUE, useFilteredControls = FALSE,
                             verbose = TRUE, loessSpan = .1, fitLoess = TRUE){  
  
  if (!is(object, "MethylSet")) stop("'object' needs to be a 'MethylSet'")
  
  # create object for methylated and unmethylated channels -------------------------
  methTmp <- getMeth(object)
  probeIDs <- rownames(methTmp)
  signals <- array(dim = c(dim(methTmp), 2))
  signals[, , 1] <- getUnmeth(object)
  signals[, , 2] <- methTmp
  frescoData <- frescoData[match(probeIDs, rownames(frescoData)), ]
  GC <- frescoData$targetGC
  
  # get set of empirical controls --------------------------------------------------
  if (useControls & !useFilteredControls){
    controls <- which(!is.na(frescoData$eControls))
    if (verbose) cat(length(controls), 'empirical control probes detected\n')
  } else if (useControls & useFilteredControls){
    controls <- which(!is.na(frescoData$eControlsFiltered))
    if (verbose) cat(length(controls), 'empirical control probes detected\n')
  }
  
  # divide probes and controls up by probe type ------------------------------------
  whichSetII <- which(frescoData$probeType == 'II')
  whichSetI <- which(frescoData$probeType == 'I') 
  
  if (useControls){ 
    whichControlsII <- intersect(whichSetII, controls)
    whichControlsI <- intersect(whichSetI, controls)
  } else {
    whichControlsII <- whichSetII
    whichControlsI <- whichSetI
  }
  
  # declare peak finding function -------------------------------------------------
  if (verbose) cat('Aligning signal intensities \n')
  
  getLowerPeak <- function(x){
    intensityDensity <- density(x)
    peaksInd <- which(diff(sign(diff(intensityDensity$y)))==-2)+1
    lowerPeak <- which.min(intensityDensity$x[peaksInd])
    return(intensityDensity$x[peaksInd[lowerPeak]])
  }
  
  # find lower peaks ---------------------------------------------------------------
  typeIpeaks <- apply(signals[whichControlsI, , ], c(2, 3), getLowerPeak)
  typeIIpeaks <- apply(signals[whichControlsII, , ], c(2, 3), getLowerPeak)
  typeIpeakMeans <- colMeans(typeIpeaks)
  typeIIpeakMeans <- colMeans(typeIIpeaks)
  
  # line up samples by their lower peaks -------------------------------------------
  signals[whichSetI, , 1] <- sweep(signals[whichSetI, , 1], 2, typeIpeaks[, 1], '-')
  signals[whichSetI, , 2] <- sweep(signals[whichSetI, , 2], 2, typeIpeaks[, 2], '-')
  signals[whichSetII, , 1] <- sweep(signals[whichSetII, , 1], 2, typeIIpeaks[, 1], '-')
  signals[whichSetII, , 2] <- sweep(signals[whichSetII, , 2], 2, typeIIpeaks[, 2], '-')
  
  # scale signals to minimize deviance from control averages -----------------------
  if (verbose) cat('Applying linear scaling factor \n')
  typeIcontrolAvg <- apply(signals[whichControlsI, , ], c(1, 3), mean, trim = .1)
  typeIIcontrolAvg <- apply(signals[whichControlsII, , ], c(1, 3), mean, trim = .1)
  
  coefsI1 <- lm(signals[whichControlsI, , 1] ~ typeIcontrolAvg[, 1] + 0)$coef
  coefsI2 <- lm(signals[whichControlsI, , 2] ~ typeIcontrolAvg[, 2] + 0)$coef
  coefsII1 <- lm(signals[whichControlsII, , 1] ~ typeIIcontrolAvg[, 1] + 0)$coef
  coefsII2 <- lm(signals[whichControlsII, , 2] ~ typeIIcontrolAvg[, 2] + 0)$coef
  
  scaledSignals <- array(dim=dim(signals))
  scaledSignals[whichSetI, , 1] <- sweep(signals[whichSetI, , 1], 2, coefsI1, '/') + typeIpeakMeans[1]
  scaledSignals[whichSetI, , 2] <- sweep(signals[whichSetI, , 2], 2, coefsI2, '/') + typeIpeakMeans[2]
  scaledSignals[whichSetII, , 1] <-sweep(signals[whichSetII, , 1], 2, coefsII1, '/') + typeIIpeakMeans[1]
  scaledSignals[whichSetII, , 2] <-sweep(signals[whichSetII, , 2], 2, coefsII2, '/') + typeIIpeakMeans[2]
  scaledSignals[scaledSignals < 0] <- 0
  
  # stop here if omitting loess fitting -------------------------------------
  if (!fitLoess){
    out <- object
    normedUnmeth <- scaledSignals[, , 1]
    normedMeth <- scaledSignalslog2NormedSignals[, , 2]
    rownames(normedUnmeth) <- rownames(normedMeth) <- probeIDs
    colnames(normedUnmeth) <- colnames(normedMeth) <- colnames(methTmp)
    
    assayDataElement(out, 'Unmeth') <- normedUnmeth
    assayDataElement(out, 'Meth') <- normedMeth
    
    out@preprocessMethod <- c(rg.norm = sprintf("fresco alignment and scaling (based on a MethylSet preprocessed as '%s'",
                                                preprocessMethod(object)[1]),
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion('IlluminaHumanMethylation450kmanifest')))
      
    return(out)
  }
  
  # compute robust experiment average ---------------------------------------------
  if (verbose) cat('Computing robust experiment-wise average\n')
  log2Centered <- log2(scaledSignals + 1)  
  sexInd <- factor(suppressWarnings(getSex(mapToGenome(object))[, 3]))
  XYind <- which(frescoData$chromosome %in% c('X', 'Y'))
  log2Standard <- apply(log2Centered, c(1, 3), mean, trim = .1)
  
  if (nlevels(sexInd) == 2){
    mInd <- which(sexInd == 'M')
    fInd <- which(sexInd == 'F')
    
    log2StandardM <- log2StandardF <- log2Standard
    log2StandardM[XYind, ] <- apply(log2Centered[XYind, mInd, ], c(1, 3), mean, trim = .1)
    log2StandardF[XYind, ] <- apply(log2Centered[XYind, fInd, ], c(1, 3), mean, trim = .1)
  }
  
  # compute deviations from average -----------------------------------------  
  if (verbose) cat('Computing deviations from average \n')
  log2Deviations <- array(dim = dim(log2Centered))
  for(kk in 1:2) 
    log2Deviations[, , kk] <- log2Centered[, , kk] - log2Standard[, kk]
  
  if (nlevels(sexInd) == 2){
    for (kk in 1:2){
      log2Deviations[XYind, mInd, kk] <- log2Centered[XYind, mInd, kk] - log2StandardM[XYind, kk]
      log2Deviations[XYind, fInd, kk] <- log2Centered[XYind, fInd, kk] - log2StandardF[XYind, kk]
    }
  }  
  
  # winsorize by probe type -------------------------------------------------
  if (useControls){
    if (verbose) cat('Winsorizing probes out of prediction range \n')
    winsorizeBySubset <- function(x, whichSet, whichControls){
      x[whichSet][which(x[whichSet] > max(x[whichControls]))] <- max(x[whichControls])
      x[whichSet][which(x[whichSet] < min(x[whichControls]))] <- min(x[whichControls])
      x[whichSet]
    }
    
    GC[whichSetI] <- winsorizeBySubset(GC, whichSetI, whichControlsI)
    GC[whichSetII] <- winsorizeBySubset(GC, whichSetII, whichControlsII)
    for (kk in 1:2){
      log2Standard[whichSetI, kk] <- winsorizeBySubset(log2Standard[, kk], whichSetI, whichControlsI)  
      log2Standard[whichSetII, kk] <- winsorizeBySubset(log2Standard[, kk], whichSetII, whichControlsII)  
    }
  }
  
  # create independent variable data frame for loess ---------------------------------
  indepVars <- data.frame(GC = GC, UMavg = log2Standard[, 1], Mavg = log2Standard[, 2])
  if(nlevels(sexInd) == 2){
    indepVarsM <- data.frame(GC = GC, UMavg = log2StandardM[, 1], Mavg = log2StandardM[, 2])
    indepVarsF <- data.frame(GC = GC, UMavg = log2StandardF[, 1], Mavg = log2StandardF[, 2])
  }
    
  # declare loess function -----------------------------------------------------------
  funLoess <- function(y, indepVars, whichControls, whichSet, smoothingParameter) {
    modelDat <- as.data.frame(cbind(y, indepVars))
    tempFit <- loess(y ~ GC * Mavg * UMavg, trace.hat = 'approx', 
                     span = smoothingParameter, 
                     modelDat, subset = whichControls)
    resids <- y[whichSet] - predict(tempFit, modelDat[whichSet, ])
    return(resids)
  }
  
  # fit loess surfaces ----------------------------------------------------------------  
  if (verbose) cat('Fitting & subtracting out loess\n')
  
  if (nlevels(sexInd) == 1){
    if (verbose) cat('Normalizing type I probes \n')
    typeInormed <- apply(log2Deviations, c(2, 3), funLoess, indepVars = indepVars, 
                         whichControls = whichControlsI, whichSet = whichSetI, smoothingParameter = loessSpan)
    
    if (verbose) cat('Normalizing type II probes \n')
    typeIInormed <- apply(log2Deviations, c(2, 3), funLoess, indepVars = indepVars, 
                          whichControls = whichControlsII, whichSet = whichSetII, smoothingParameter = loessSpan)
    
    # recombine into log2NormedDevs
    log2NormedDevs <- array(dim = dim(log2Deviations))
    log2NormedDevs[whichSetI, , ] <- typeInormed
    log2NormedDevs[whichSetII, , ] <- typeIInormed  
  }
  
  if (nlevels(sexInd) == 2){
    if (verbose) cat('Normalizing type I probes \n')
    # type I
    typeInormedM <- apply(log2Deviations[, mInd, ], c(2, 3), funLoess, indepVars = indepVarsM, 
                          whichControls = whichControlsI, whichSet = whichSetI, smoothingParameter = loessSpan)
    
    typeInormedF <- apply(log2Deviations[, fInd, ], c(2, 3), funLoess, indepVars = indepVarsF, 
                          whichControls = whichControlsI, whichSet = whichSetI, smoothingParameter = loessSpan)
    
    # type II
    if (verbose) cat('Normalizing type II probes \n')
    typeIInormedM <- apply(log2Deviations[, mInd, ], c(2, 3), funLoess, indepVars = indepVarsM, 
                           whichControls = whichControlsII, whichSet = whichSetII, smoothingParameter = loessSpan)
    
    cat('Normalizing type II probes \n')
    typeIInormedF <- apply(log2Deviations[, fInd, ], c(2, 3), funLoess, indepVars = indepVarsF, 
                           whichControls = whichControlsII, whichSet = whichSetII, smoothingParameter = loessSpan)
    
    # recombine into log2NormedDevs
    log2NormedDevs <- array(dim = dim(log2Deviations))
    log2NormedDevs[whichSetI, mInd, ] <- typeInormedM
    log2NormedDevs[whichSetII, mInd, ] <- typeIInormedM
    log2NormedDevs[whichSetI, fInd, ] <- typeInormedF
    log2NormedDevs[whichSetII, fInd, ] <- typeIInormedF
    
  }
  
  # compute normalized log2 signals --------------------------------------------------
  log2NormedSignals <- array(dim = dim(log2Centered))  
  if (nlevels(sexInd) == 1)
    for (kk in 1:2) log2NormedSignals[, , kk] <- log2NormedDevs[, , kk] + log2Standard[, kk]
  
  if (nlevels(sexInd) == 2){
    for(kk in 1:2){
      log2NormedSignals[, mInd, kk] <- log2NormedDevs[, mInd, kk] + log2StandardM[, kk]
      log2NormedSignals[, fInd, kk] <- log2NormedDevs[, fInd, kk] + log2StandardF[, kk]
    }
  }
  
  # create new MethylSet for output ------------------------------------------------
  out <- object
  normedUnmeth <- 2^log2NormedSignals[, , 1]
  normedMeth <- 2^log2NormedSignals[, , 2]
  rownames(normedUnmeth) <- rownames(normedMeth) <- probeIDs
  colnames(normedUnmeth) <- colnames(normedMeth) <- colnames(methTmp)
  
  assayDataElement(out, 'Unmeth') <- normedUnmeth
  assayDataElement(out, 'Meth') <- normedMeth
  
  out@preprocessMethod <- c(rg.norm = sprintf("fresco alignment, scaling, and surface fitting (based on a MethylSet preprocessed as '%s'",
                                            preprocessMethod(object)[1]),
                            minfi = as.character(packageVersion('minfi')),
                            manifest = as.character(packageVersion('IlluminaHumanMethylation450kmanifest')))
  
  out
}

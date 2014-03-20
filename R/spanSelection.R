#' Get fit statistics R^2, AICC, and GCV for loess fits
#' 
#' @param object \code{MethylSet} object
#' @param useControls Should empirical controls be used to align and fit loess surfaces?
#' @param loessSpan Supply span for fitting loess surface
#' @param sdThreshold Threshold to filter empirical controls by standard deviation

returnFitStats <-function(object, useControls = TRUE, loessSpan = .15, 
                          sdThreshold = .15, verbose = FALSE){  
  
  if (loessSpan > 1 | loessSpan < 0) stop("loessSpan must be between zero and one")
  
  # create object for methylated and unmethylated channels -------------------------
  methTmp <- getMeth(object)
  probeIDs <- rownames(methTmp)
  signals <- array(dim = c(dim(methTmp), 2))
  signals[, , 1] <- getUnmeth(object)
  signals[, , 2] <- methTmp
  frescoData <- frescoData[match(probeIDs, rownames(frescoData)), ]
  GC <- frescoData$targetGC
  
  # get set of empirical controls --------------------------------------------------
  if (useControls){
    probeSD <- apply(getBeta(object), 1, sd)
    controls <- which(!is.na(frescoData$eControls) & probeSD < sdThreshold)
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
  
  # find lower peaks ---------------------------------------------------------------
  if (verbose) cat('Aligning signal intensities \n')
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
  
  scaledSignals <- array(dim = dim(signals))
  scaledSignals[whichSetI, , 1] <- sweep(signals[whichSetI, , 1], 2, coefsI1, '/') + typeIpeakMeans[1]
  scaledSignals[whichSetI, , 2] <- sweep(signals[whichSetI, , 2], 2, coefsI2, '/') + typeIpeakMeans[2]
  scaledSignals[whichSetII, , 1] <-sweep(signals[whichSetII, , 1], 2, coefsII1, '/') + typeIIpeakMeans[1]
  scaledSignals[whichSetII, , 2] <-sweep(signals[whichSetII, , 2], 2, coefsII2, '/') + typeIIpeakMeans[2]
  scaledSignals[scaledSignals < 0] <- 0
  
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
  
  # winsorize by probe type ----------------------------------------------------------
  if (useControls){
    if (verbose) cat('Winsorizing probes out of prediction range \n')
    
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
  
  # fit loess surfaces ----------------------------------------------------------------  
  if (verbose) cat('Fitting & subtracting out loess\n')
  
  if (nlevels(sexInd) == 1){
    if (verbose) cat('Normalizing type I probes \n')
    typeInormed <- apply(log2Deviations, c(2, 3), funLoessSS, 
                         indepVars = indepVars, whichControls = whichControlsI,
                         whichSet = whichSetI, smoothingParameter = loessSpan)
    
    if (verbose) cat('Normalizing type II probes \n')
    typeIInormed <- apply(log2Deviations, c(2, 3), funLoessSS, 
                          indepVars = indepVars, whichControls = whichControlsII, 
                          whichSet = whichSetII, smoothingParameter = loessSpan)
    
    estimatedErrorVar <- list(typeInormed, typeIInormed)
    return(estimatedErrorVar)
  }
  
  if (nlevels(sexInd) == 2){
    if (verbose) cat('Normalizing type I probes \n')
    # type I
    typeInormedM <- apply(log2Deviations[, mInd, ], c(2, 3), funLoessSS,
                          indepVars = indepVarsM, whichControls = whichControlsI, 
                          whichSet = whichSetI, smoothingParameter = loessSpan)
    
    typeInormedF <- apply(log2Deviations[, fInd, ], c(2, 3), funLoessSS, 
                          indepVars = indepVarsF, whichControls = whichControlsI, 
                          whichSet = whichSetI, smoothingParameter = loessSpan)
    
    # type II
    if (verbose) cat('Normalizing type II probes \n')
    typeIInormedM <- apply(log2Deviations[, mInd, ], c(2, 3), funLoessSS, 
                           indepVars = indepVarsM, whichControls = whichControlsII, 
                           whichSet = whichSetII, smoothingParameter = loessSpan)
    
    typeIInormedF <- apply(log2Deviations[, fInd, ], c(2, 3), funLoessSS, 
                           indepVars = indepVarsF, whichControls = whichControlsII, 
                           whichSet = whichSetII, smoothingParameter = loessSpan)
    
    typeInormed <- typeIInormed <- array(dim = c(3, dim(log2Deviations)[2], 2))
    typeInormed[, mInd, ] <- typeInormedM
    typeInormed[, fInd, ] <- typeInormedF
    typeIInormed[, mInd, ] <- typeIInormedM
    typeIInormed[, fInd, ] <- typeIInormedF
    
    estimatedErrorVar <- list(typeInormed, typeIInormed)
    return(estimatedErrorVar)
  }
  
}


# declare loess function 
funLoessSS <- function(y, indepVars, whichControls, whichSet, smoothingParameter) {
  modelDat <- as.data.frame(cbind(y, indepVars))
  tempFit <- loess(y ~ GC * Mavg * UMavg, trace.hat = 'approx', span = smoothingParameter, 
                   modelDat, subset = whichControls)
  traceL <- tempFit$trace.hat
  sigma2 <- sum(tempFit$residuals^2)/(tempFit$n - 1)
  aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(tempFit$n - traceL - 2)
  gcv <- tempFit$n * sigma2/(tempFit$n - traceL)^2    
  r2 <- cor(tempFit$y, tempFit$fitted)^2
  return(c(aicc, gcv, r2))
}

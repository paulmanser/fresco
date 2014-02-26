#' Get fit statistics R^2, AICC, and GCV for loess fits
#' 
#' @param object a \code{MethylSet} object
#' @param useControls Should empirical controls be used to align and fit loess surfaces?
#' @param loessSpan Supply vector of possible spans for fitting loess surface
#' @param sdThreshold Threshold to filter empirical controls by standard deviation
#' 
#' @export plotFitStats


plotFitStats <- function(object, useControls = TRUE, 
                         loessSpan = seq(.05, .95, .15), sdThreshold = .1){
  
  if (!is(object, "MethylSet")) stop("'object' needs to be a 'MethylSet'")
  
  data(frescoData)
  
  fitstats <- list()
  
  # generate fit statistics -------------------------------------------------
  for(ii in 1:length(loessSpan)){
    fitstats[[ii]] <- returnFitStats(object, useControls = useControls, 
                                     loessSpan = loessSpan[ii], sdThreshold = sdThreshold)
    cat(ii, 'of', length(sParm), '\n\n')
  }
  
  pdf('fresco-fit-stats.pdf')
  par(mfrow = c(2,2))
  
  # generate gcv curves -----------------------------------------------------
  chType <- c('UM', 'M')
  statType <- c('AICC', 'GCV', 'R^2')
  for(thisStat in 1:3){
    for(probeType in 1:2){
      for(channelType in 1:2){
      
      gcvCurves <- NULL
      for(ii in 1:length(fitstats))
        gcvCurves <- cbind(gcvCurves, fitstats[[ii]][[probeType]][2, , channelType])
      
      matplot(loessSpan, t(gcvCurves), type='l', 
              main = paste('Type ', probeType, ' : ',
                           chType[channelType], ' channel : ',
                           statType[thisStat], sep = '')
              )
      }
    }
  }  
  dev.off()
}
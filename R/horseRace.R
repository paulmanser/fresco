
horseRace <- function(object, batchVarName = NULL, 
                      covariateNames = NULL, covariateTypes = NULL){
  
  if (!is(object, "RGChannelSet"))
    stop("object needs to be a 'RGChannelSet'")
  
  if (is.null(batchVarName) & is.null(covariateNames))
    stop("Please provide variable name for either batch or covariates of interest")
  
  if (!is.null(covariateNames) & is.null(covariateTypes))
    stop("Please provide corresponding vector of covariate types: either 
         'categorical' or 'continuous'")
  
  if (!is.null(covariateTypes) & any(!covariateTypes %in% c('categorical', 'continuous')))
    stop("covariateTypes must be one of 'categorical' or 'continuous'")
  
  # normalize ---------------------------------------------------------
  normList <- list()
  normList$rawMSet  <- preprocessIllumina(updateObject(object))
  normList$fresco15 <- preprocessFresco(normList$rawMSet, loessSpan = .15, sdThreshold = .15)
  normList$fresco50 <- preprocessFresco(normList$rawMSet, loessSpan = .5, sdThreshold = .15)
  normList$fresco85 <- preprocessFresco(normList$rawMSet, loessSpan = .85, sdThreshold = .15)
  normList$frescoNL <- preprocessFresco(normList$rawMSet, fitLoess = FALSE, sdThreshold = .15)
  normList$quantile <- preprocessQuantile(mapToGenome(object))
  normList$funnorm  <- preprocessFunnorm(object)

  # test for batch effects ---------------------------------------------
  if (!is.null(batchVarName)){
    f.results <- lapply(normList, .batchTest)
    roc.results <- lapply(f.results, function(x) .rocComp(na.omit(x[, 2])))
  }
  
  # p-value ecdf
  plot(0, 0, xlim = c(0, 1), ylim = c(0, 1),
       xlab = 'P-value', ylab = 'ECDF',
       main = 'P-value ECDF for Batch Effects')
  abline(0, 1, lty = 3)
  
  for(ii in 1:length(roc.results))
    lines(seq(0, 1, .01), roc.results[[ii]], col = ii)
  legend('bottomright', legend = names(roc.results), fill = 1:length(roc.results))
  
  # look at power for covariates of interest --------------------------------
  
  if (!is.null(covariateNames)){
    for(ii in 1:length(covariateNames)){
      
      if (is(covariateNames[ii], 'categorical')){
        f.results <- lapply(normList, .catTest)
        roc.results <- lapply(f.results, function(x) .rocComp(na.omit(x[, 2]), BH.adj = TRUE))
      }
      
      if (is(covariateNames[ii], 'continuous')){
        f.results <- lapply(normList, .contTest)
        roc.results <- lapply(f.results, function(X) .rocComp(na.omit(x[, 2]), BH.adj = TRUE))
      }
      
    # plot
    plot(0, 0, xlim = c(0, 1), ylim = c(0, 1),
         xlab = 'FDR Threshold', ylab = 'Prop sig at given FDR',
         main = paste('Power for', covariateNames[iii]))
    abline(0, 1, lty = 3)
    
    for(ii in 1:length(roc.results))
      lines(seq(0, 1, .01), roc.results[[ii]], col = ii)
    legend('bottomright', legend = names(roc.results), fill = 1:length(roc.results))
    }    
  }
  
}


.rocComp <- function(x, eval = seq(0, 1, .01), BH.adj = FALSE){
  if (BH.adj) x <- p.adjust(x, method = 'BH')
  sapply(eval, function(z) sum(x < z) / length(x))
}

.batchTest <- function(x) rowFtests(getBeta(x), factor(pData(x)[, batchVarName]))
.catTest <- function(x) rowFtests(getBeta(x), factor(pData(x)[, factor(covariateNames[ii])]))
.contTest <- function(x){
  cont.cov <- pData(x)[, covariateNames[ii]]
  apply(getBeta(x), 1, function(z) biglm(z ~ cont.cov))
}





horseRace <- function(object, batchVarName = NULL, 
                      covariateNames = NULL, covariateTypes = NULL,
                      compositeF=FALSE){
  
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
  normList$Raw       <- preprocessIllumina(updateObject(object))
  normList$FRESCO_15 <- preprocessFresco(normList$Raw, loessSpan = .15, sdThreshold = .1)
  normList$FRESCO_50 <- preprocessFresco(normList$Raw, loessSpan = .5, sdThreshold = .1)
  normList$FRESCO_85 <- preprocessFresco(normList$Raw, loessSpan = .85, sdThreshold = .1)
  normList$FRESCO_NL <- preprocessFresco(normList$Raw, fitLoess = FALSE, sdThreshold = .1)
  normList$SQN       <- preprocessQuantile(mapToGenome(object))
  normList$Funnorm   <- preprocessFunnorm(object)
  normList$Noob      <- preprocessNoob(object)

  # test for batch effects ---------------------------------------------
  if (!is.null(batchVarName)){
    f.results <- lapply(normList, .catTest, cvn=batchVarName)
    roc.results <- lapply(f.results, function(x) .rocComp(na.omit(x[, 2])))
    
    # p-value ecdf
    plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type='n',
         xlab = 'P-value', ylab = 'ECDF',
         main = 'P-value ECDF for Batch Effects')
    abline(0, 1, lty = 3)
    
    for(ii in 1:length(roc.results))
      lines(seq(0, 1, .01), roc.results[[ii]], col = ii)
    legend('bottomright', legend = names(roc.results), fill = 1:length(roc.results))
  }
  
  # look at power for covariates of interest --------------------------------
  
  if (!is.null(covariateNames) & !compositeF){
    for(ii in 1:length(covariateNames)){
      
      if (covariateTypes[ii] == 'categorical'){
        f.results <- lapply(normList, .catTest, cvn=covariateNames[ii])
        roc.results <- lapply(f.results, function(x) .rocComp(na.omit(x[, 2]), BH.adj = TRUE))
      }
      
      if (covariateTypes[ii] == 'continuous'){
        f.results <- lapply(normList, .contTest, cvn=covariateNames[ii])
        roc.results <- lapply(f.results, function(X) .rocComp(na.omit(x[, 2]), BH.adj = TRUE))
      }
      
      # plot
      plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type='n',
           xlab = 'FDR Threshold', ylab = 'Prop sig at given FDR',
           main = paste('Significant Differences for', covariateNames[ii]))
      abline(0, 1, lty = 3)
      
      for(ii in 1:length(roc.results))
        lines(seq(0, 1, .01), roc.results[[ii]], col = ii)
      legend('bottomright', legend = names(roc.results), fill = 1:length(roc.results))
    }    
  }
  
  # look at composite F-scores for covariates of interest --------------------
  
  if (!is.null(covariateNames) & compositeF){
    for(ii in 1:length(covariateNames)){
      
      if (covariateTypes[ii] == 'categorical'){
        
        # compute anova sum of squares and get degrees of freedom
        f.results <- lapply(normList, .compSS, cvn=covariateNames[ii])
        p <- nlevels(factor(pData(normList[[1]])[, covariateNames[ii]]))
        n <- ncol(normList[[1]])
        
        # re-order f-results because mapToGenome re-orders the CpGs
        for(jj in 2:length(f.results)){
          f.results[[jj]] <- f.results[[jj]][match(rownames(normList[[1]]), rownames(normList[[jj]])) ,]
        }
        
        p.vals <- list()        
        
        for(jj in 2:length(normList)){
          # compute composite F scores and p-values
          comp.pvals <- comp.f.stats <- matrix(nr=nrow(normList[[1]]), nc=3)
          colnames(comp.pvals) <- colnames(comp.f.stats) <- c('orig', 'ssr_raw', 'sse_raw')
          comp.f.stats[, 1] <- (f.results[[1 ]][, 1] / (p-1)  ) / (f.results[[1 ]][, 2] / (n-p)  )
          comp.f.stats[, 2] <- (f.results[[1 ]][, 1] / (p-1)  ) / (f.results[[jj]][, 2] / (n-p)  )
          comp.f.stats[, 3] <- (f.results[[jj]][, 1] / (p-1)  ) / (f.results[[1 ]][, 2] / (n-p)  )
          comp.pvals <- pf(comp.f.stats, df1=p-1, df2=n-p, lower.tail=FALSE)
          p.vals[[(jj-1)]] <- comp.pvals
        }
        names(p.vals) <- names(normList)[-1]
        
        # plot
        axis.lims <- -log10(unlist(p.vals))
        axis.lims <- max(axis.lims[which(is.finite(axis.lims))])
        
        par(mfcol = c(2, 5), mar=c(5, 5, 5, 3))

        for(jj in c(1, 4:7)){
          plot( -log10(p.vals[[jj]][, 1]), -log10(p.vals[[jj]][, 2]),
               pch=16, cex=.2, col=rgb(0,0,1,alpha=.4),
               xlab = 'Original F Statistic', ylab = expression('F'['Err']),
               main=names(p.vals)[jj], xlim=c(0, axis.lims), ylim=c(0, axis.lims),
               cex.axis = 1.7, cex.lab = 1.7, cex.main=1.7)
          abline(0,1,col='red')
          
          plot( -log10(p.vals[[jj]][, 1]),  -log10(p.vals[[jj]][, 3]),
               pch=16, cex=.2, col=rgb(0,0,1,alpha=.4),
               xlab = 'Original F Statistic', ylab = expression('F'['ES']),
               xlim=c(0, axis.lims), ylim=c(0, axis.lims),
               cex.axis = 1.7, cex.lab = 1.7, cex.main=1.7)
          abline(0,1,col='red')
          
        }
      }
      
      if (covariateTypes[ii] == 'continuous'){
        cat('Not yet supported')
      }
      
    }  
  }
  
  
}


.rocComp <- function(x, eval = seq(0, 1, .01), BH.adj = FALSE){
  if (BH.adj) x <- p.adjust(x, method = 'BH')
  sapply(eval, function(z) sum(x < z) / length(x))
}

.catTest <- function(x, cvn) rowFtests(getBeta(x), factor(pData(x)[, cvn]))
.contTest <- function(x, cvn){
  cont.cov <- as.numeric(pData(x)[, cvn])
  apply(getBeta(x), 1, function(z) biglm(z ~ cont.cov))
}

.compSS <- function(x, cvn){
  cat.cov <- factor(pData(x)[, cvn])
  betas <- getBeta(x)
  fac.means <- matrix(nr=nrow(betas), nc=nlevels(cat.cov))
  res.mat <- matrix(nr=nrow(betas), nc=ncol(betas))
  mean.vec <- rowMeans(betas)
  
  for(ii in 1:nlevels(cat.cov)){
    fac.level <- which(cat.cov == levels(cat.cov)[ii])
    fac.means[, ii] <- rowMeans(betas[, fac.level])
    res.mat[, fac.level] <- betas[, fac.level] - fac.means[, ii]
  }

  sse <- rowSums(res.mat^2)
  ssr <- rowSums(sweep((fac.means - mean.vec)^2, 2, table(cat.cov), '*'))
  
  out <- cbind(ssr, sse)
  out
}


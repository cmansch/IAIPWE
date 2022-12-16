

getBoundariesChi <- function(filedir, alpha=0.05){
  # initialize the Z value for boundaries at quantile for univariate testing
  bounds <- list()
  for (mod in c('aa', 'cd', 'ipw')){
    # get the z correlation
    corMatrix <- as.matrix(read.csv(paste0(filedir, "/Zcorrelation", mod, ".csv"), row.names=1))

    # CHANGE JANUARY 05 2022 TO ACCOUNT FOR IOTA INFORMATION PROPORTION
    varEsts <- rlist::list.load(paste0(filedir, "/varEstimateResults.Rdata"))
    varname <- paste("aipw", mod, sep="_")
    if (mod == "ipw"){
      varname <- mod
    }
    # n1/n2
    # propint <- mean(unlist(lapply(varEsts, function(x) x$t1[[varname]]$n / x$t2[[varname]]$n )))
    # varVector <- as.matrix(read.csv(paste0(filedir, "/estimatedvariance", mod, ".csv"), row.names=1))
    # v1 <- varVector[1:(length(varVector)/2)] 
    # v2 <- varVector[(length(varVector)/2 + 1):length(varVector)] 
    # iota <- sqrt((mean(v1)/mean(v2))  * (1/propint))
    i2 <- mean(unlist(lapply(varEsts, function(x) (mean(diag(x$t1[[varname]]$variance)) / x$t1[[varname]]$n)  / 
                               (mean(diag(x$t2[[varname]]$variance)) / x$t2[[varname]]$n) ) ) )
    # END CHANGE JANUARY 5 2022
    
    fixedBound <- findChiSquareBound(corr=corMatrix[(dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1]),
                                                   (dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])], 
                                    alpha=0.05, mu=NULL, B=10000, tol=1e-4, 
                                    startvalue=1, lambda=0.1, seed=1, K=1,
                                    iota=c(1))
    pocockBound <- findChiSquareBound(corr=corMatrix, alpha=0.05, mu=NULL, B=10000, tol=1e-4, 
                                  startvalue=1, lambda=0.1, seed=1, K=2,
                                  iota=c(1,1))
    
    ofBound <- findChiSquareBound(corr=corMatrix, alpha=0.05, mu=NULL, B=10000, tol=1e-4, 
                                   startvalue=1, lambda=0.1, seed=1, K=2,
                                   iota=c(i2,1))
    
    bounds[[mod]] <- list('Fixed' = fixedBound,
                          'Pocock'  = pocockBound,
                          'OF' = ofBound,
                          'iota' = i2)
  }
  out_location <-  paste0(filedir, "/stoppingBoundaries.Rdata")
  rlist::list.save(bounds, file = out_location)
  return()
}
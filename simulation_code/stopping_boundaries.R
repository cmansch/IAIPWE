
getBoundaries <- function(filedir, alpha=0.05){
  # initialize the Z value for boundaries at quantile for univariate testing
  initZvalue <- qnorm(1-alpha)
  bounds <- list()
  for (mod in c('aa', 'cd', 'ipw')){
    # get the z correlation
    corMatrix <- as.matrix(read.csv(paste0(filedir, "/Zcorrelation", mod, ".csv"), row.names=1))
    # find the Pocock type boundaries (i.e. calpha, calpha)
    pocockBound <- findPocockBoundary(corr=corMatrix, alpha=alpha, startvalue = initZvalue)
    # find OBrien Fleming type boundaries (e.g. sqrt(2)calpha, calpha)
    
    # CHANGE JANUARY 05 2022 TO ACCOUNT FOR IOTA INFORMATION PROPORTION
    varEsts <- rlist::list.load(paste0(filedir, "/varEstimateResults.Rdata"))
    varname <- paste("aipw", mod, sep="_")
    if (mod == "ipw"){
      varname <- mod
    }
    # n1/n2
    # get all iotas
    L <- length(varEsts[[1]][["t1"]][[varname]]$value)
    S <- length(grep("c_t", names(varEsts[[1]])))   # number of analyses
    varVector <- as.matrix(read.csv(paste0(filedir, "/estimatedvariance", mod, ".csv"), row.names=1))
    varSmap <- rep(1:S, each=L)
    iota <- c()
    for (s in 1:(S-1)){
      propint <-  mean(unlist(lapply(varEsts, function(x) (x[[paste0("t",s)]][[varname]]$n + x[[paste0("c_t",s)]]$n) / 
                                       (x[[paste0("t",S)]][[varname]]$n  + x[[paste0("c_t",S)]]$n    )  )))
      vs <- varVector[varSmap == s]  # interim 
      vS <- varVector[varSmap == S]  # final
      iota <- c(iota, sqrt( (mean(vs) / mean(vS)) * (1/propint) ))
    }
        
    # propint <- mean(unlist(lapply(varEsts, function(x) x$t1[[varname]]$n / x$t2[[varname]]$n )))
    # varVector <- as.matrix(read.csv(paste0(filedir, "/estimatedvariance", mod, ".csv"), row.names=1))
    # v1 <- varVector[1:(length(varVector)/2)] 
    # v2 <- varVector[(length(varVector)/2 + 1):length(varVector)] 
    # iota <- sqrt((mean(v1)/mean(v2))  * (1/propint))
    # END CHANGE JANUARY 5 2022
    
    ofBound <- findOFBoundary(corr=corMatrix, alpha=alpha, 
                              startvalue = initZvalue, iota=iota, L=L)
    # find the boundary that would be needed if the analysis is only done
    # at the final time point
    fixedBound <- findPocockBoundary(corr=corMatrix[(dim(corMatrix)[1] - L +1) :(dim(corMatrix)[1]),
                                                    (dim(corMatrix)[1] - L +1) :(dim(corMatrix)[1])],
                                     alpha=alpha, startvalue = initZvalue)
    bounds[[mod]] <- list('Fixed' = fixedBound,
                          'Pocock'  = pocockBound,
                          'OF' = ofBound,
                          'iota' = iota)
  }
  out_location <-  paste0(filedir, "/stoppingBoundaries.Rdata")
  rlist::list.save(bounds, file = out_location)
  return()
}
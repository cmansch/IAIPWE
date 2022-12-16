
# assumes that there are 2 analyses. 

#### Find the sample sizes required for alternative true values ####
# HaltDelta is the differences in the regimes from null under alternative.
# for example, if regime 1 and 4 are 2 higher, it would be c(2,0,0,2).
# lambda is the sample size step size initially
# N0 is the sample size to start at
# propint is the proportion of inidividuals through the trial when the interim
#  analysis is to be performed; depricated
# propint is now estimated by what generated the variance results, because that 
#  ratio is what the standard errors are valid for (i.e. when the analysis for 
#  interim time is completed)
getSampleSizesChi <- function(filedir, power, HaltDelta=HaltDelta, lambda=10, N0=50, propint=0.5){
  # read in bounds
  allbounds <- rlist::list.load(paste0(filedir, "/stoppingBoundaries.Rdata"))
  varEsts <- rlist::list.load(paste0(filedir, "/varEstimateResults.Rdata"))
  
  samples <- list()
  
  for (mod in c('aa', 'cd', 'ipw')){
    # find the proper interim analysis alternatives based on the average sample 
    # size at the interim analysis
    varname <- paste("aipw", mod, sep="_")
    if (mod == "ipw"){
      varname <- mod
    }
    # n1/n2
    propint <- mean(unlist(lapply(varEsts, function(x) x$t1[[varname]]$n / x$t2[[varname]]$n )))
    
    bounds <- allbounds[[mod]]
    # get the z correlation
    corMatrix <- as.matrix(read.csv(paste0(filedir, "/Zcorrelation", mod, ".csv"), row.names=1))
    # get the estimated variance
    varVector <- as.matrix(read.csv(paste0(filedir, "/estimatedvariance", mod, ".csv"), row.names=1))
    # print(varVector)
    # print(corMatrix)
    # print(HaltDelta)
    # print(power)
    # print(bounds$Pocock$bound)
    # find the Pocock type sample size (i.e. calpha, calpha)
    
    ## get the sample sizes; single, pocock, of
    fixedN <- findChiSquareSampleFixed(corr=corMatrix[((dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])),
                                                     ((dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1]))], 
                                      pow=power, 
                                      cov=varVector[(dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])], 
                                      vp=HaltDelta, calpha=bounds$Fixed$bound,
                                      lambda=10, N0=50, seed=1, 
                                      K=1, B=10000)
    
    pocockN <- findChiSquareSample(corMatrix, pow=power, cov=varVector, vp=as.vector(rep(HaltDelta, 2)), 
                               calpha=bounds$Pocock$bound, 
                               interim=propint,
                               lambda=lambda, N0=N0, iota=c(1,1), seed=1, 
                               K=2, B=10000) 
    
    ofN <- findChiSquareSample(corMatrix, pow=power, cov=varVector, vp=as.vector(rep(HaltDelta, 2)), 
                                calpha=bounds$OF$bound, interim=propint,
                                lambda=lambda, N0=N0, iota=c(bounds$iota,1), seed=1, 
                                K=2, B=10000)
    
    samples[[mod]] <- list('Fixed' = fixedN,
                           'Pocock'  = pocockN,
                           'OF' = ofN,
                           'iota' = bounds$iota)
  }
  out_location <-  paste0(filedir, "/sampleSizes.Rdata")
  rlist::list.save(samples, file = out_location)
  print("The sample sizes required for the alternative proposed are: ")
  print(samples)
  return()
}
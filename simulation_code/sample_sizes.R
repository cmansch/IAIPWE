

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
getSampleSizes <- function(filedir, power, HaltDelta=HaltDelta, lambda=10, N0=50, 
                           propint=0.5, S=2){  # S=2 to give option of multiple 
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
    L <- length(varEsts[[1]][["t1"]][[varname]]$value)
    S <- length(grep("c_t", names(varEsts[[1]])))  # number of analyses
    # n1/n2
    # propint <- mean(unlist(lapply(varEsts, function(x) x$t1[[varname]]$n / x$t2[[varname]]$n )))
    propint <- c()
    for (s in 1:(S-1)){
      temp <-  mean(unlist(lapply(varEsts, function(x) (x[[paste0("t",s)]][[varname]]$n + x[[paste0("c_t",s)]]$n) / 
                                    (x[[paste0("t",S)]][[varname]]$n  + x[[paste0("c_t",S)]]$n    )  )))
      propint <- c(propint, temp)
    }
    
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
    pocockN <- findPocockSample(corr=corMatrix, pow=power, cov=varVector, 
                                vp=as.vector(rep(HaltDelta, S)),
                                bound=bounds$Pocock$bound, interim=propint)
    
    # CHANGE MADE ON JANUARY 5 2022 TO ACCOUNT FOR INFORMATION PROPORTION OF OF BOUNDARY
    # get the estimated proportion of information 
    # v1 <- varVector[1:(length(varVector)/2)] 
    # v2 <- varVector[(length(varVector)/2 + 1):length(varVector)] 
    # iota <- sqrt((mean(v1)/mean(v2))  * (1/propint))
    
    # find OBrien Fleming type sample size (e.g. sqrt(2)calpha, calpha)
    # use bounds$iota since this is now available following stopping boundaries
    ofN <- findOFSample(corMatrix, power, varVector, vp=as.vector(rep(HaltDelta, S)),
                        bound=bounds$OF$bound, interim=propint, iota=bounds$iota) 
    
    # find the sample size that would be needed if the analysis is only done
    # at the final time point
    fixedN <- findFixedSample(corMatrix[(dim(corMatrix)[1] - L +1) :(dim(corMatrix)[1]),
                                        (dim(corMatrix)[1] - L +1) :(dim(corMatrix)[1])],
                              pow=power,
                              cov=varVector[(dim(corMatrix)[1] - L +1) :(dim(corMatrix)[1])],
                              vp=HaltDelta, bound=bounds$Fixed$bound,
                              lambda=lambda, N0=N0)
    
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
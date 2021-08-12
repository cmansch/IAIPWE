## simulation to compare ipw, aa-aipw, cd-aipw
## we use design 1 data to do this.

if(substr(getwd(), 2, 5) == 'home'){
  rm(list=ls())
}


## source the files containing any necessary functions
## and libraries
library(dplyr)
library(foreach)
library(doParallel)
library(smartim)

## define the setting
Bvar <- 2000  # simulations to estimate the variance
nvar <- 1000  # sample size for variance estimation
s2 <- 100  # variance of outcome Y
regimes <- list(c(0,0), c(0,1), c(1,0), c(1,1))  # list of embedded regimes in data
steps <- 2  # number of decisions K
delta <- 10  # the value under the null hypothesis
vpvar <- 1
vpsim <- 2
Bsim <- 2000
H0 <- c(10,10,10,10)
ncores <- 20


#### simulate the variance ####

estimateAllMethod <- function(n, loopno, vp, regimes, steps, s2){
  # make data
  dat <- design1(vp=vp, n=n, seed=n*loopno, s2=s2)

  # get analysis times
  t1 <- median(dat$t3)
  t2 <- max(dat$t3)

  # give the q functions
  q2 <- modelObj::buildModelObj(model = ~ a1 + a2 + a1:a2,
                                solver.method = 'lm',
                                predict.method = 'predict.lm')
  q1 <- modelObj::buildModelObj(model = ~ a1,
                                solver.method = 'lm',
                                predict.method = 'predict.lm')
  q_list <- list(q1, q2)
  # get list of regimes to be used in final and all available
  L <- length(regimes)
  regime_all <- regimelist(embRegimes = regimes,
                           dat = dat)

  ### estimate results for AIPW AA
  im_mid_aa <- vAIPW(dat=dat, ti=t1,
                     steps=steps,
                     q_list=q_list,
                     regime_list=regime_all,
                     feasibleSetsIndicator=NULL)
  im_final <- vAIPW(dat=dat, ti=t2,
                    steps=steps,
                    q_list=q_list,
                    regime_list=regime_all,
                    feasibleSetsIndicator=NULL)
  # subset the data for completed data and ipw approach.
  # To do completed data only, we subset the dat data set. Otherwise the
  # VAIPW code would use partial information

  dat_mid <- dat[dat$t3 <= t1,]
  regime_mid <- regimelist(embRegimes = regimes,
                           dat = dat_mid)
  ### AIPW CD
  im_mid_cd <- vAIPW(dat=dat_mid, ti=t1,
                     steps=steps,
                     q_list=q_list,
                     regime_list=regime_mid,
                     feasibleSetsIndicator=NULL)

  ### now do IPW estimator
  im_mid_ipw <- vIPW(dat_mid, regime_mid)
  im_fin_ipw <- vIPW(dat, regime_all)

  #### Create the output for each of the three methods ####

  ## AIPW AA
  # interim
  im_mid_aa_val <- unlist(im_mid_aa$Values)
  im_mid_aa_dimCov <- dim(im_mid_aa$Covariance)[1]
  im_mid_aa_varOfVal <- im_mid_aa$Covariance[(im_mid_aa_dimCov-L+1):im_mid_aa_dimCov,
                                             (im_mid_aa_dimCov-L+1):im_mid_aa_dimCov ]
  im_mid_aa_se <- sqrt(diag(im_mid_aa_varOfVal) / im_mid_aa$ns)
  # final
  im_final_val <- unlist(im_final$Values)
  im_final_dimCov <- dim(im_final$Covariance)[1]
  im_final_varOfVal <- im_final$Covariance[(im_final_dimCov-L+1):im_final_dimCov,
                                           (im_final_dimCov-L+1):im_final_dimCov ]
  im_final_se <- sqrt(diag(im_final_varOfVal) / im_final$ns)

  # output
  aipw_aa <- list('t1' = list('value' = im_mid_aa_val,
                              'variance' = im_mid_aa_varOfVal,
                              'se' = im_mid_aa_se,
                              'seed' = n*loopno,
                              'n' = im_mid_aa$ns),
                  't2' = list('value' = im_final_val,
                              'variance' = im_final_varOfVal,
                              'se' = im_final_se,
                              'seed' = n*loopno,
                              'n' = n))
  ## AIPW CD
  # interim
  im_mid_cd_val <- unlist(im_mid_cd$Values)
  im_mid_cd_dimCov <- dim(im_mid_cd$Covariance)[1]
  im_mid_cd_varOfVal <- im_mid_cd$Covariance[(im_mid_cd_dimCov-L+1):im_mid_cd_dimCov,
                                             (im_mid_cd_dimCov-L+1):im_mid_cd_dimCov ]
  im_mid_cd_se <- sqrt(diag(im_mid_cd_varOfVal) / im_mid_cd$ns)

  aipw_cd <- list('t1' = list('value' = im_mid_cd_val,
                              'variance' = im_mid_cd_varOfVal,
                              'se' = im_mid_cd_se,
                              'seed' = n*loopno,
                              'n' = im_mid_cd$ns),
                  't2' = list('value' = im_final_val,
                              'variance' = im_final_varOfVal,
                              'se' = im_final_se,
                              'seed' = n*loopno,
                              'n' = n))
  ## IPW
  im_mid_ipw_val <- unlist(im_mid_ipw$Values)
  im_mid_ipw_dimCov <- dim(im_mid_ipw$Covariance)[1]
  im_mid_ipw_varOfVal <- im_mid_ipw$Covariance[(im_mid_ipw_dimCov-L+1):im_mid_ipw_dimCov,
                                               (im_mid_ipw_dimCov-L+1):im_mid_ipw_dimCov ]
  im_mid_ipw_se <- sqrt(diag(im_mid_ipw_varOfVal) / im_mid_ipw$ns)
  # final
  im_fin_ipw_val <- unlist(im_fin_ipw$Values)
  im_fin_ipw_dimCov <- dim(im_fin_ipw$Covariance)[1]
  im_fin_ipw_varOfVal <- im_fin_ipw$Covariance[(im_fin_ipw_dimCov-L+1):im_fin_ipw_dimCov,
                                               (im_fin_ipw_dimCov-L+1):im_fin_ipw_dimCov ]
  im_fin_ipw_se <- sqrt(diag(im_fin_ipw_varOfVal) / im_fin_ipw$ns)

  # output
  ipw <- list('t1' = list('value' = im_mid_ipw_val,
                          'variance' = im_mid_ipw_varOfVal,
                          'se' = im_mid_ipw_se,
                          'seed' = n*loopno,
                          'n' = im_mid_ipw$ns),
              't2' = list('value' = im_fin_ipw_val,
                          'variance' = im_fin_ipw_varOfVal,
                          'se' = im_fin_ipw_se,
                          'seed' = n*loopno,
                          'n' = n))
  return(list('aipw_aa' = aipw_aa,
              'aipw_cd' = aipw_cd,
              'ipw' = ipw))

}

estimateVariance <- function(n, loops, vp, regimes, steps, s2,
                             delta, out_location=NULL, cores=20){
  print('Simulating data')
  print(Sys.time())
  # results <- list()
  # for (i in 1:loops){
  #   results[[i]] <- estimateAllMethod(n=n, loopno=i, vp=vp, regimes=regimes, steps=steps, s2=s2)
  # }
  registerDoParallel(cores=cores)

  results <- foreach(n=rep(n, loops), loopno=c(1:loops), vp=rep(vp,loops),
                     regimes=replicate(n=loops, expr=regimes, simplify = FALSE),
                     steps=rep(steps, loops), s2=rep(s2, loops)) %dopar% {
                       estimateAllMethod(n, loopno, vp, regimes, steps, s2)
                     }
  rlist::list.save(results, file = paste0(out_location, "/varEstimateResults.Rdata"))

  print('Estimating the variance')
  print(Sys.time())

  # list input should have structure [[seed]][[t1/t2]][[results from seed/time]]
  singleVar <- function(list, delta, out_location, mod){
    z1s <- lapply(list, function(x) (x$t1$value - delta)/x$t1$se)
    z2s <- lapply(list, function(x) (x$t2$value - delta)/x$t2$se)
    # find the covariance of the z1s, z2s
    zall <- cbind(matrix(unlist(z1s), byrow=TRUE, nrow=length(z1s)),
                  matrix(unlist(z2s), byrow=TRUE, nrow=length(z2s)))
    # and empirical variances
    varvalues <- colMeans(matrix(unlist(lapply(list,
                                             function(x) c(diag(x$t1$variance),
                                                           diag(x$t2$variance)))),
                               ncol=ncol(zall), byrow=TRUE))
    # save output
    write.csv(as.data.frame(cov(zall)), paste0(out_location, "/Zcovariance", mod, ".csv"))
    write.csv(as.data.frame(cor(zall)), paste0(out_location, "/Zcorrelation", mod, ".csv"))
    write.csv(as.data.frame(varvalues), paste0(out_location, "/estimatedvariance", mod,".csv"))
    return()
  }
  ### now find all single variances.
  ## aa
  singleVar(lapply(results, function(x) x$aipw_aa), delta, out_location, 'aa')
  ## cd
  singleVar(lapply(results, function(x) x$aipw_cd), delta, out_location, 'cd')
  ## ipw
  singleVar(lapply(results, function(x) x$ipw), delta, out_location, 'ipw')
  return()
}


#### find stopping boundaries given variance ####

## currently only works for SMARTs with 2 interim analyses
getBoundaries <- function(filedir, alpha=0.05){
  # initialize the Z value for boundaries at quantile for univariate testing
  initZvalue <- qnorm(1-alpha)
  bounds <- list()
  for (mod in c('aa', 'cd', 'ipw')){
    # get the z correlation
    corMatrix <- as.matrix(read.csv(paste0(filedir, "/Zcorrelation", mod, ".csv"), row.names=1))
    # find the Pocock type boundaries (i.e. calpha, calpha)
    pocockBound <- findPocockBoundary(corr=corMatrix, alpha=alpha, startvalue = initZvalue)
    # find OBrien Fleming type boundaries (i.e. sqrt(2)calpha, calpha)
    ofBound <- findOFBoundary(corr=corMatrix, alpha=alpha, startvalue = initZvalue)
    # find the boundary that would be needed if the analysis is only done 
    # at the final time point
    fixedBound <- findPocockBoundary(corr=corMatrix[(dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1]),
                                               (dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])],
                                     alpha=alpha, startvalue = initZvalue)
    bounds[[mod]] <- list('Fixed' = fixedBound,
                          'Pocock'  = pocockBound,
                          'OF' = ofBound)
  }
  out_location <-  paste0(filedir, "/stoppingBoundaries.Rdata")
  rlist::list.save(bounds, file = out_location)
  return()
}

# HaltDelta is the differences in the regimes from null under alternative. 
# for example, if regime 1 and 4 are 2 higher, it would be c(2,0,0,2). 
# lambda is the sample size step size initially
# N0 is the sample size to start at
# propint is the proportion of inidividuals through the trial when the interim 
# analysis is to be performed
getSampleSizes <- function(filedir, power, HaltDelta, lambda=10, N0=50, propint=0.5){
  # read in bounds
  allbounds <- rlist::list.load(paste0(filedir, "/stoppingBoundaries.Rdata"))
  
  samples <- list()
  
  for (mod in c('aa', 'cd', 'ipw')){
    bounds <- allbounds[[mod]]
    # get the z correlation
    corMatrix <- as.matrix(read.csv(paste0(filedir, "/Zcorrelation", mod, ".csv"), row.names=1))
    # get the estimated variance 
    varVector <- as.matrix(read.csv(paste0(filedir, "/estimatedvariance", mod, ".csv"), row.names=1))
    
    # find the Pocock type sample size (i.e. calpha, calpha)
    pocockN <- findPocockSample(corMatrix, power, varVector, vp=as.vector(rep(HaltDelta, 2)),
                     bound=bounds$Pocock$bound, interim=propint)
    
    # find OBrien Fleming type sample size (i.e. sqrt(2)calpha, calpha)
    ofN <- findOFSample(corMatrix, power, varVector, vp=rep(HaltDelta, 2),
                 bound=bounds$Pocock$bound, interim=propint)
    
    # find the sample size that would be needed if the analysis is only done 
    # at the final time point
    fixedN <- findFixedSample(corMatrix[((dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])), 
                                        ((dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1]))], 
                                 pow=power,
                                 cov=varVector[(dim(corMatrix)[1]/2+1) :(dim(corMatrix)[1])],
                                 vp=HaltDelta, bound=bounds$Fixed$bound,
                                 lambda=lambda, N0=N0)
    
    samples[[mod]] <- list('Fixed' = fixedN,
                          'Pocock'  = pocockN,
                          'OF' = ofN)
  }
  out_location <-  paste0(filedir, "/sampleSizes.Rdata")
  rlist::list.save(samples, file = out_location)
  
  return()
}

#### Simulation for average sample size and simulated power ####
# allN is the list of sample sizes for power
# these should be [[aa, cd, ipw]][[Pocock, OF, Fixed]][[N, power]]

getK1ResultsList <- function(list, seed){
  values <- unlist(list$Values)
  L <- length(values)
  dimCov <- dim(list$Covariance)[1]
  varOfValues <- list$Covariance[(dimCov-L+1):dimCov,
                                               (dimCov-L+1):dimCov ]
  ses <- sqrt(diag(varOfValues) / list$ns)
  
  # output
  return( list('value' = values,
               'variance' = varOfValues,
               'se' = ses,
               'seed' = seed,
               'n' = list$ns)) 
}

# getK2ResultsList <- function(list1, list2, seed){
#   t1res <- getK1ResultsList(list1, seed)
#   t2res <- getK1ResultsList(list2, seed)
#   return(list('t1' = t1res,
#               't2' = t2res))
# }

ipwSingle <- function(dat, regimes, seed, bound, truen){
  # get regimes
  regime_all <- regimelist(embRegimes = regimes, dat = dat)
  # run it
  ipw_fixed <- vIPW(dat, regime_all)
  # get results list
  ipw_fixed_res <- getK1ResultsList(ipw_fixed, seed)
  # check if rejected
  # get if any were rejected
  # then find the sum of how many were rejected
  # and if it is greater than 0, then the reject is 1
  ipw_fixed_rej <- (sum((((ipw_fixed_res$value - H0) / ipw_fixed_res$se) > 
                           bound) > 0) > 0)*1
  # output effective sample size (min if early rejected, max else)
  # ipw_fixed_sample_size <- ipw_fixed_res$n
  
  return(list('results' = ipw_fixed_res,
              'rejectInd' = ipw_fixed_rej,
              'sample' = truen))
}

simulationOnce <- function(maxN, allN, bounds, regimes, vp, loopno, s2, H0,
                           steps){
  ## generate the maximum data set
  # make data
  seed <- maxN*loopno
  dat <- design1(vp=vp, n=maxN, seed=seed, s2=s2)
  # give the q functions
  q2 <- modelObj::buildModelObj(model = ~ a1 + a2 + a1:a2,
                                solver.method = 'lm',
                                predict.method = 'predict.lm')
  q1 <- modelObj::buildModelObj(model = ~ a1,
                                solver.method = 'lm',
                                predict.method = 'predict.lm')
  q_list <- list(q1, q2)
  # get list of regimes to be used in final and all available
  L <- length(regimes)
  
  ### ipw ########################################################
  ipwres <- list()
  ipwN <- allN[['ipw']]
  ipwBounds <- bounds$ipw
  ## Fixed ##
  ipw_dat_f <- dat[1:ipwN$Fixed$N,]
  ipwres[['Fixed']] <- ipwSingle(ipw_dat_f, regimes, seed, 
                                 ipwBounds$Fixed$bound, nrow(ipw_dat_f))
  
  ## Pocock ##
  ipw_dat_p <- dat[1:ipwN$Pocock$N[[2]],]
  t1 <- median(ipw_dat_p$t3)
  ipw_dat_p_mid <- ipw_dat_p[ipw_dat_p$t3 <= t1,]
  nmid <- nrow(ipw_dat_p[ipw_dat_p$t1 <= t1,])  # how many individuals are there when we analyze t1

  ipwres[['Pocock_Mid']] <- ipwSingle(ipw_dat_p_mid, regimes, seed, 
                                      ipwBounds$Pocock$bound, truen=nmid)
  ipwres[['Pocock_Fin']] <- ipwSingle(ipw_dat_p, regimes, seed, 
                                      ipwBounds$Pocock$bound, truen=nrow(ipw_dat_p))
  # did IPW reject early?
  # did it reject at least once?
  # what is the sample size for this study?
  ipwres[['Pocock']] <- list('earlyReject' =  ipwres[['Pocock_Mid']]$rejectInd,
                             'anyReject' =  ((ipwres[['Pocock_Mid']]$rejectInd + 
                                                ipwres[['Pocock_Fin']]$rejectInd) > 0)*1,
                             'sampleSize' = (ipwres[['Pocock_Mid']]$rejectInd * ipwres[['Pocock_Mid']]$sample) + 
                               ((1-ipwres[['Pocock_Mid']]$rejectInd)* ipwres[['Pocock_Fin']]$sample)
  )
  
  ## OBrien Fleming ##
  ipw_dat_of <- dat[1:ipwN$OF$N[[2]],]
  t1 <- median(ipw_dat_of$t3)
  ipw_dat_of_mid <- ipw_dat_of[ipw_dat_of$t3 <= t1,]
  nmid <- nrow(ipw_dat_of[ipw_dat_of$t1 <= t1,]) 
  
  ipwres[['OF_Mid']] <- ipwSingle(ipw_dat_of_mid, regimes, seed, 
                                      ipwBounds$OF$bound*sqrt(2), truen=nmid)
  ipwres[['OF_Fin']] <- ipwSingle(ipw_dat_of, regimes, seed, 
                                      ipwBounds$OF$bound, truen=nrow(ipw_dat_of))
  # did IPW reject early?
  # did it reject at least once?
  # what is the sample size for this study?
  ipwres[['OF']] <- list('earlyReject' =  ipwres[['OF_Mid']]$rejectInd,
                             'anyReject' =  ((ipwres[['OF_Mid']]$rejectInd + 
                                                ipwres[['OF_Fin']]$rejectInd) > 0)*1,
                             'sampleSize' = (ipwres[['OF_Mid']]$rejectInd * ipwres[['OF_Mid']]$sample) + 
                               ((1-ipwres[['OF_Mid']]$rejectInd)* ipwres[['OF_Fin']]$sample)
  )
  
  ### AIPW AA #################################################################
  aares <- list()
  aaN <- allN[['aa']]
  aaBounds <- bounds$aa
  ### Fixed
  ## subset data
  aa_dat <- dat[1:aaN$Fixed$N,]
  # get times
  t1 <- median(aa_dat$t3)
  t2 <- max(aa_dat$t3)
  # get regimes
  regime_all <- regimelist(embRegimes = regimes,
                           dat = aa_dat)
  # run
  aa_fin <- vAIPW(dat=aa_dat, ti=t2,
                    steps=steps,
                    q_list=q_list,
                    regime_list=regime_all,
                    feasibleSetsIndicator=NULL)
  aa_fixed_res <- getK1ResultsList(aa_fin, seed)
  aa_fixed_rej <- (sum((((aa_fixed_res$value - H0) / aa_fixed_res$se) > 
                           aaBounds$Fixed$bound) > 0) > 0)*1
  # output effective sample size (min if early rejected, max else)
  aares[['Fixed']] <- list('results' = aa_fixed_res,
                           'rejectInd' = aa_fixed_rej,
                           'sample' = aa_fixed_res$n)
  
  ### Pocock
  ## subset data
  aa_dat <- dat[1:aaN$Pocock$N[[2]],]
  # get times
  t1 <- median(aa_dat$t3)
  t2 <- max(aa_dat$t3)
  # get regimes
  regime_all <- regimelist(embRegimes = regimes,
                           dat = aa_dat)
  # run
  aa_mid <- vAIPW(dat=aa_dat, ti=t1,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  aares[['Pocock_Mid']] <- getK1ResultsList(aa_mid, seed)
  aa_fin <- vAIPW(dat=aa_dat, ti=t2,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  aares[['Pocock_Fin']] <- getK1ResultsList(aa_fin, seed)
  # check for early reject
  aa_p_mid <- (sum((((aares[['Pocock_Mid']]$value - H0) / aares[['Pocock_Mid']]$se) > 
                        aaBounds$Pocock$bound) > 0) > 0)*1
  aa_p_fin <- (sum((((aares[['Pocock_Fin']]$value - H0) / aares[['Pocock_Fin']]$se) > 
                        aaBounds$Pocock$bound) > 0) > 0)*1
  # give results
  aares[['Pocock']] <- list('earlyReject' = aa_p_mid ,
                            'anyReject'= ((aa_p_mid + aa_p_fin)>0)*1,
                            'sampleSize'= (aa_p_mid*aares[['Pocock_Mid']]$n) + 
                              ((1-aa_p_mid)*aares[['Pocock_Fin']]$n ) )
  
  ### OBrien Fleming 
  aa_dat <- dat[1:aaN$OF$N[[2]],]
  # get times
  t1 <- median(aa_dat$t3)
  t2 <- max(aa_dat$t3)
  # get regimes
  regime_all <- regimelist(embRegimes = regimes,
                           dat = aa_dat)
  # run
  aa_mid <- vAIPW(dat=aa_dat, ti=t1,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  aares[['OF_Mid']] <- getK1ResultsList(aa_mid, seed)
  aa_fin <- vAIPW(dat=aa_dat, ti=t2,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  aares[['OF_Fin']] <- getK1ResultsList(aa_fin, seed)
  # check for early reject
  aa_p_mid <- (sum((((aares[['OF_Mid']]$value - H0) / aares[['OF_Mid']]$se) > 
                      aaBounds$OF$bound*sqrt(2)) > 0) > 0)*1
  aa_p_fin <- (sum((((aares[['OF_Fin']]$value - H0) / aares[['OF_Fin']]$se) > 
                      aaBounds$OF$bound) > 0) > 0)*1
  # give results
  aares[['OF']] <- list('earlyReject' = aa_p_mid ,
                            'anyReject'= ((aa_p_mid + aa_p_fin)>0)*1,
                            'sampleSize'= (aa_p_mid*aares[['OF_Mid']]$n) + 
                              ((1-aa_p_mid)*aares[['OF_Fin']]$n ) )
  
  ### AIPW CD #################################################################
  cdres <- list()
  cdN <- allN[['cd']]
  cdBounds <- bounds$cd
  ### Fixed
  ## subset data
  cd_dat <- dat[1:cdN$Fixed$N,]
  # get times
  t1 <- median(cd_dat$t3)
  t2 <- max(cd_dat$t3)
  # get regimes
  regime_all <- regimelist(embRegimes = regimes,
                           dat = cd_dat)
  # run
  cd_fin <- vAIPW(dat=cd_dat, ti=t2,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  cd_fixed_res <- getK1ResultsList(cd_fin, seed)
  cd_fixed_rej <- (sum((((cd_fixed_res$value - H0) / cd_fixed_res$se) > 
                          cdBounds$Fixed$bound) > 0) > 0)*1
  # output effective sample size (min if early rejected, max else)
  cdres[['Fixed']] <- list('results' = cd_fixed_res,
                           'rejectInd' = cd_fixed_rej,
                           'sample' = cd_fixed_res$n)
  
  ### Pocock
  # similar to aa but now we actually subset the data at interim time
  ## subset data
  cd_dat <- dat[1:cdN$Pocock$N[[2]],]
  # get times
  t1 <- median(cd_dat$t3)
  t2 <- max(cd_dat$t3)
  # subset mid data
  cd_dat_mid <- cd_dat[cd_dat$t3 <= t1,]
  nmid <- nrow(cd_dat[cd_dat$t1 <= t1,]) 
  # get regimes
  regime_mid <- regimelist(embRegimes = regimes,
                           dat = cd_dat_mid)
  regime_all <- regimelist(embRegimes = regimes,
                           dat = cd_dat)
  # run
  cd_mid <- vAIPW(dat=cd_dat_mid, ti=t1,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_mid,
                  feasibleSetsIndicator=NULL)
  cdres[['Pocock_Mid']] <- getK1ResultsList(cd_mid, seed)
  cd_fin <- vAIPW(dat=cd_dat, ti=t2,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  cdres[['Pocock_Fin']] <- getK1ResultsList(cd_fin, seed)
  # check for early reject
  cd_p_mid <- (sum((((cdres[['Pocock_Mid']]$value - H0) / cdres[['Pocock_Mid']]$se) > 
                      cdBounds$Pocock$bound) > 0) > 0)*1
  cd_p_fin <- (sum((((aares[['Pocock_Fin']]$value - H0) / cdres[['Pocock_Fin']]$se) > 
                      cdBounds$Pocock$bound) > 0) > 0)*1
  # give results
  cdres[['Pocock']] <- list('earlyReject' = cd_p_mid ,
                            'anyReject'= ((cd_p_mid + cd_p_fin)>0)*1,
                            'sampleSize'= (cd_p_mid*nmid) + 
                              ((1-cd_p_mid)*nrow(cd_dat) ) )
  ### OBrien Fleming
  cd_dat <- dat[1:cdN$OF$N[[2]],]
  # get times
  t1 <- median(cd_dat$t3)
  t2 <- max(cd_dat$t3)
  # subset mid data
  cd_dat_mid <- cd_dat[cd_dat$t3 <= t1,]
  nmid <- nrow(cd_dat[cd_dat$t1 <= t1,]) 
  # get regimes
  regime_mid <- regimelist(embRegimes = regimes,
                           dat = cd_dat_mid)
  regime_all <- regimelist(embRegimes = regimes,
                           dat = cd_dat)
  
  # run
  # run
  cd_mid <- vAIPW(dat=cd_dat_mid, ti=t1,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_mid,
                  feasibleSetsIndicator=NULL)
  cdres[['OF_Mid']] <- getK1ResultsList(cd_mid, seed)
  cd_fin <- vAIPW(dat=cd_dat, ti=t2,
                  steps=steps,
                  q_list=q_list,
                  regime_list=regime_all,
                  feasibleSetsIndicator=NULL)
  cdres[['OF_Fin']] <- getK1ResultsList(cd_fin, seed)
  # check for early reject
  cd_p_mid <- (sum((((cdres[['OF_Mid']]$value - H0) / cdres[['OF_Mid']]$se) > 
                      cdBounds$OF$bound*sqrt(2)) > 0) > 0)*1
  cd_p_fin <- (sum((((aares[['OF_Fin']]$value - H0) / cdres[['OF_Fin']]$se) > 
                      cdBounds$OF$bound) > 0) > 0)*1
  # give results
  cdres[['OF']] <- list('earlyReject' = cd_p_mid ,
                        'anyReject'= ((cd_p_mid + cd_p_fin)>0)*1,
                        'sampleSize'= (cd_p_mid*nmid) + 
                          ((1-cd_p_mid)*nrow(cd_dat) ) )
  
  ## return save values, variances, and zscores
  
  return(list('ipw' = ipwres,
              'aa' = aares,
              'cd' = cdres))
  
}


simulationParallel <- function(filedir, regimes, vp, s2, H0, loops, steps, cores){
  sampleSizes <- rlist::list.load(paste0(filedir, "/sampleSizes.Rdata"))
  maxSampleSize <- max(unlist(lapply(sampleSizes, function(x) c(x[[1]]$N, x[[2]]$N, x[[3]]$N))))
  bounds <- rlist::list.load(paste0(filedir, "/stoppingBoundaries.Rdata"))
  ## run simulations in parallel
  
  registerDoParallel(cores=cores)
  
  results <- foreach(maxN=rep(maxSampleSize, loops), 
                     allN=replicate(n=loops, expr=sampleSizes, simplify=FALSE),
                     bounds=replicate(n=loops, expr=bounds, simplify=FALSE),
                     regimes=replicate(n=loops, expr=regimes, simplify = FALSE),
                     vp=rep(vp,loops),
                     loopno=c(1:loops), 
                     s2=rep(s2, loops),
                     H0=replicate(n=loops, expr=H0, simplify = FALSE),
                     steps=rep(steps, loops)) %dopar% {
       simulationOnce(maxN, allN, bounds, regimes, vp, 
                      loopno, s2, H0, steps)
  }
  rlist::list.save(results, file = paste0(filedir, "/simulationData.Rdata"))
  
  return()
}

summaryStats <- function(filedir){
  print('This will calculate the power and average sample size')
  results <- rlist::list.load(file = paste0(filedir, "/simulationData.Rdata"))
  bounds <- rlist::list.load(file = paste0(filedir, "/stoppingBoundaries.Rdata"))
  # find the average sample size
  aaN <- c(mean(unlist(lapply(results, function(x) x$aa$Fixed$sample))),
           mean(unlist(lapply(results, function(x) x$aa$Pocock$sampleSize))),
           mean(unlist(lapply(results, function(x) x$aa$OF$sampleSize))) )         
  cdN <- c(mean(unlist(lapply(results, function(x) x$cd$Fixed$sample))),
           mean(unlist(lapply(results, function(x) x$cd$Pocock$sampleSize))),
           mean(unlist(lapply(results, function(x) x$cd$OF$sampleSize))) )         
  ipwN<- c(mean(unlist(lapply(results, function(x) x$ipw$Fixed$sample))),
           mean(unlist(lapply(results, function(x) x$ipw$Pocock$sampleSize))),
           mean(unlist(lapply(results, function(x) x$ipw$OF$sampleSize))) )         
  
  # find the proportion of early rejects 
  aaEarly <- c(mean(unlist(lapply(results, function(x) x$aa$Fixed$rejectInd))),
               mean(unlist(lapply(results, function(x) x$aa$Pocock$earlyReject))),
               mean(unlist(lapply(results, function(x) x$aa$OF$earlyReject))) )         
  cdEarly <- c(mean(unlist(lapply(results, function(x) x$cd$Fixed$rejectInd))),
               mean(unlist(lapply(results, function(x) x$cd$Pocock$earlyReject))),
               mean(unlist(lapply(results, function(x) x$cd$OF$earlyReject))) )         
  ipwEarly <- c(mean(unlist(lapply(results, function(x) x$ipw$Fixed$rejectInd))),
                mean(unlist(lapply(results, function(x) x$ipw$Pocock$earlyReject))),
                mean(unlist(lapply(results, function(x) x$ipw$OF$earlyReject))) )         
  # find the proportion of total rejects
  aaAny <- c(mean(unlist(lapply(results, function(x) x$aa$Fixed$rejectInd))),
               mean(unlist(lapply(results, function(x) x$aa$Pocock$anyReject))),
               mean(unlist(lapply(results, function(x) x$aa$OF$anyReject))) )         
  cdAny <- c(mean(unlist(lapply(results, function(x) x$cd$Fixed$rejectInd))),
               mean(unlist(lapply(results, function(x) x$cd$Pocock$anyReject))),
               mean(unlist(lapply(results, function(x) x$cd$OF$anyReject))) )         
  ipwAny <- c(mean(unlist(lapply(results, function(x) x$ipw$Fixed$rejectInd))),
                mean(unlist(lapply(results, function(x) x$ipw$Pocock$anyReject))),
                mean(unlist(lapply(results, function(x) x$ipw$OF$anyReject))) ) 
  
  summaryStat <- list('aaN' = aaN, 
                      'cdN' = cdN, 
                      'ipwN' = ipwN, 
                      'aaEarly' = aaEarly, 
                      'cdEarly' = cdEarly, 
                      'ipwEarly'= ipwEarly, 
                      'aaAny' = aaAny, 
                      'cdAny' = cdAny, 
                      'ipwAny' = ipwAny)
  
  # check the values
  # apply((matrix(unlist(lapply(results, function(x) x$aa$Fixed$results$value)), ncol=4, byrow=TRUE)), 2, mean)
  # apply((matrix(unlist(lapply(results, function(x) (x$aa$Pocock_Mid$value - 10)/ x$aa$Pocock_Mid$se)), ncol=4, byrow=TRUE)), 2, mean)
  # apply((matrix(unlist(lapply(results, function(x) x$aa$Pocock_Fin$value)), ncol=4, byrow=TRUE)), 2, mean)
  # 
  
  print(summaryStat)
  
  rlist::list.save(summaryStat, file = paste0(filedir, "/simulationSummary.Rdata"))
  
  return()
}



#### main function ####
# prevdate should be formatted as "2021-08-10"
main <- function(varest=TRUE, boundest=TRUE, sampleest=TRUE, simyn=TRUE,
                 unpacksim=TRUE, newdir=TRUE, prevdate=NULL){
  # handle output director
  if(newdir){
    toppath <- paste0(getwd(), '/sim_output/', Sys.Date())
    if(!dir.exists(toppath)){
      dir.create(toppath)
      sendmailR::sendmail(from = '<design1>',
                          to = '<cmansch@ncsu.edu>',
                          subject = 'Server Update',
                          msg = paste('Directory made successfully') )
    }
  } else if(!is.null(prevdate)){
    toppath <- paste0(getwd(), '/sim_output/', prevdate)
    if (!dir.exists(toppath)){
      print('Directory given does not exist')
      break
    }
  } else{
    print('No new directory created and previous date is missing.')
    break
  }
  
  # estimate the variance of the data under each of the estimators
  if(varest){
    estimateVariance(n=nvar, loops=Bvar, vp=vpvar, regimes=regimes,
                     steps=steps, s2=s2, delta=delta, out_location=toppath,
                     cores=ncores)
  }
  # get stopping boundaries under different Value Patterns
  if (boundest){
    print('Getting boundaries')
    print(Sys.time())
    getBoundaries(toppath, alpha=0.05)
  }
  # get sample sizes under stopping boundaries
  if (sampleest){
    print('Getting sample sizes')
    print(Sys.time())
    getSampleSizes(toppath, power=0.8, HaltDelta=c(0,0,0,2), lambda=10, N0=50, 
                   propint=0.5)
  }
  # simulate for sample sizes and return the observed power and average sample size
  if (simyn){
    print('Starting simulation')
    print(Sys.time())
    simulationParallel(filedir=toppath, regimes=regimes, vp=vpsim, 
                       s2=s2, H0=H0, loops=Bsim, steps=steps, cores=ncores)
    print('Simulations over')
    print(Sys.time())
  }
  if (unpacksim){
    print('Summarizing results')
    print(Sys.time())
    summaryStats(filedir=toppath)
  }
  
  sendmailR::sendmail(from = '<design1>',
                      to = '<cmansch@ncsu.edu>',
                      subject = 'Server Update',
                      msg = paste('Code has finished running. Check the server for results') )

}

#### run main ####
# run the main function if on the server
if(substr(getwd(), 2, 5) == 'home'){
  main(varest=TRUE, sampleest=TRUE, boundest=TRUE, simyn=TRUE,
       unpacksim = TRUE, newdir=TRUE, prevdate=NULL)  # "2021-08-12"
  #2021-08-10
}


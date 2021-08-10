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
source(paste(getwd(), "value.R" , sep="/"))
source(paste(getwd(), "vipw.R" , sep="/"))
source(paste(getwd(), "design1data.R" , sep="/"))
source(paste(getwd(), "regimes.R" , sep="/"))
source(paste(getwd(), "findPocockOB.R" , sep="/"))

## define the setting
Bvar <- 2000  # simulations to estimate the variance
nvar <- 1000  # sample size for variance estimation
s2 <- 100  # variance of outcome Y
regimes <- list(c(0,0), c(0,1), c(1,0), c(1,1))  # list of embedded regimes in data
steps <- 2  # number of decisions K
delta <- 10  # the value under the null hypothesis
vpvar <- 1

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
                               ncol=nrow(zall), byrow=TRUE))
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


#### main function ####
main <- function(varest=TRUE, saveyn=TRUE, newdir=TRUE){
  if(newdir){
    toppath <- paste0(getwd(), '/sim_output/', Sys.Date())
    if(!dir.exists(toppath)){
      dir.create(toppath)
      sendmailR::sendmail(from = '<design1>',
                          to = '<cmansch@ncsu.edu>',
                          subject = 'Server Update',
                          msg = paste('Directory made successfully') )
    }
  }


  # estimate the variance of the data under each of the estimators
  if(varest){
    estimateVariance(n=nvar, loops=Bvar, vp=vpvar, regimes=regimes,
                     steps=steps, s2=s2, delta=delta, out_location=toppath,
                     cores=20)
  }
  # get stopping boundaries under different Value Patterns

  # get sample sizes under stopping boundaries

  # simulate for sample sizes and return the observed power and average sample size

  sendmailR::sendmail(from = '<design1>',
                      to = '<cmansch@ncsu.edu>',
                      subject = 'Server Update',
                      msg = paste('Code has finished running. Check the server for results') )

}

#### run main ####
# run the main function if on the server
if(substr(getwd(), 2, 5) == 'home'){
  main()
}


## find the stopping boundaries and sample sizes


# this will actually run 5 different estimates, but it should be faster than it is

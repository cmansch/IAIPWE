#### Make data function ####
# still 47.5
# r2p depricated
# expects a 2 character number for vp
# the first is the value pattern
# the second characterizes the proportion of individuals enrolled at time t=500
design3 <- function(vp, n, seed=1, s2=100, a1p=0.5, a2p=0.5, r2p = 0.5){
  set.seed(seed)
  vp <- as.character(vp)
  time_pattern <- substr(vp, 2, 2)
  vp <- substr(vp, 1, 1)
  # get the parameter values associated with the value pattern given by user
  if (vp==1){
    # to vary x higher, use 0.02, 3, 3
    # vpbeta <- c(10, 0, 0, 0, 0.5, 12.5, 12.5, 0, 0)  # all equal
    vpbeta <- c(5, 0, 0, 0, 0.5, -10, 0.4, 0, 0)  # all equal
  } else if (vp==2) {
    vpbeta <- c(8, -3, 0, 0, 0.5, -10, 0.4, 0, 0)  # a1==0 better uniformly
  } else if (vp==3) {
    vpbeta <- c(10, -5, 0, 0, 0.5, -10, 0.4, 0, 0)   #a1==0 even better
  } else if (vp==4) {
    vpbeta <- c(7, -2, 0, 0, 0.5, -10, 0.4, 0, 0)   #a1==0 even better
  }
  
  ## test vps
  a1 <- c(0,0,0,0,1,1,1,1)
  r2 <- c(0,0,1,1,0,0,1,1)
  a2 <- c(0,1,0,1,0,1,0,1)
  x11 <- 47.5
  x12 <- 0.5
  x21 <- 47.5 * 1.25
  vs <- cbind(1, a1, a2, a1*a2, x11, x12, x21, a1*x11, a1*x12) %*% vpbeta
  mean(vs[c(1,3)])
  mean(vs[c(1,4)])
  mean(vs[c(2,3)])
  mean(vs[c(2,4)])
  mean(vs[c(5,7)])
  mean(vs[c(5,8)])
  mean(vs[c(6,7)])
  mean(vs[c(6,8)])
  
  # randomly assign treatments to each individual
  a1 <- rbinom(n=n, size=1, prob=a1p)
  # create the error for individuals
  ei <- rnorm(n=n, mean=0, sd=sqrt(s2))
  # create the observed outcome based
  # aDesignMatrix <- cbind(1, a1, a2, a1*a2)
  
  # now the design specific for covariates (mean 0)
  x11 <- rnorm(n=n, mean=47.5, sd=8) # mean 47.5 
  x12 <- rbinom(n=n, size=1, p=0.5) # mean 0.5 # gender
  x21 <- rnorm(n=n, mean=x11 * 1.25, sd=3)
  # interim evaluated yi, using vpbeta[2] for the coefficient of a1
  # doing so correlates x21 with a1 and covariates
  # we comment this out, but use for the second simulation
  # x21 <- 10 + 0.02*x11 - 3*x12 + x13 + vpbeta[2]*a1
  # for covariate to be mean zero, we use c(0.02, -3, 1), do this to show
  # equivalence with previous approaches
  
  # r2 <- rbinom(n=n, size=1, prob=r2p)
  expit <- function(x){ exp(x) / (1 + exp(x))}
  # want the prob to have approximately mean 0
  # summary(r2_prob)
  r2_prob <- expit(0.01*x11 + 0.02*x12 - 0.008*x21)
  # plot(density(r2_prob)) # target 50% outcomes 
  r2 <- rbinom(n=n, size=1, prob=r2_prob)
  
  a2 <- rbinom(n=n, size=1, prob=a2p)  
  
  # y <- (aDesignMatrix %*% vpbeta) + cbind(x11, x12, x21) %*% c(0.02, 3, 3) + a1*x11 + a1*x12 + ei
  y <- (cbind(1, a1, a2, a1*a2, x11, x12, x21, a1*x11, a1*x12) %*% vpbeta) + ei
  
  # create the times of arrival for patients
  # these are done using the fixed times given to control for variance due
  # to enrollment process.
  # We multiple the process by 1000 to have the unit be 'days'.
  # time patterns are based on proportion to stage k at time 500. 
  # pattern \ K | 1 | 2 | 3 |
  #   1         | 10| 10| 30|
  #   2         | 10| 20| 20|
  #   3         | 10| 30| 10|
  #   4         | 20| 10| 20|
  #   5         | 20| 20| 10|
  #   6         | 30| 10| 10|
  # time patterns are based on proportion to stage k at time 700.   # used in study
  # pattern \ K | 1 | 2 | 3 |
  #   1         | 10| 10| 50|
  #   2         | 10| 20| 40|
  #   3         | 10| 30| 30|
  #   4         | 20| 10| 40|
  #   5         | 20| 20| 30|
  #   6         | 30| 10| 30|
  # proportion completed, then at stage 2, then enrolled
  if (time_pattern==1){
    psplit <- c(0, 0.5, 0.6, 0.7, 1)  # 1, 1, 5
  }
  if (time_pattern==2){
    psplit <- c(0, 0.4, 0.6, 0.7, 1)  # 1, 2, 4
  }
  if (time_pattern==3){
    psplit <- c(0, 0.3, 0.6, 0.7, 1)  # 1, 3, 3
  }
  if (time_pattern==4){
    psplit <- c(0, 0.4, 0.5, 0.7, 1)  # 2, 1, 4
  }
  if (time_pattern==5){
    psplit <- c(0, 0.3, 0.5, 0.7, 1)  # 2, 2, 3
  }
  if (time_pattern==6){
    psplit <- c(0, 0.3, 0.4, 0.7, 1)  # 3, 1, 3
  }
  
  u <- runif(n=n, min=0, max=1)
  buckets <- cut(u, breaks=psplit, labels=1:4)
  
  t1 <- (buckets == 1) * runif(n=n, min=0, max=500) + 
    (buckets == 2) * runif(n=n, min=501, max=600) + 
    (buckets == 3) * runif(n=n, min=601, max=700) + 
    (buckets == 4) * runif(n=n, min=701, max=1000) 
  
  t2 <- t1 + 0.1*1000
  t3 <- t2 + 0.1*1000
  
  # if (time_pattern==1){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  # # the jump is done using the 0.1*1000 where .1 jumps 10%, .2 jumps 20%, etc.
  # if (time_pattern==2){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.2*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  # if (time_pattern==3){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.3*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  # if (time_pattern==4){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.2*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  # if (time_pattern==5){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.2*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.2*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  # if (time_pattern==6){
  #   t1 <- round(runif(n, 0, 1), 3)*1000
  #   t2 <- t1 + 0.3*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  #   t3 <- t2 + 0.1*1000 # + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # }
  
  dat <- as.data.frame(cbind(a1, a2, y, r2, x11, x12, x21, t1, t2, t3))
  names(dat)[3] <- "y"  
  return(dat)
}


#### Regime Indicator function ####
regimelist3 <- function(embRegimes, dat){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]])
  
  for(r in c(1:nRegimes)){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)
    
    regimen <- embRegimes[[r]]
    for (k in 1:K){
      if (paste0("r", k) %in% colnames(dat)) {
        # then we have trt for non respond followed by trt for resp
        regime[,k] <- regimen[[k]][1]*(1-dat[,paste0("r", k)]) + 
          regimen[[k]][2]*(dat[,paste0("r", k)])
        
      } else{
        regime[,k] <- regimen[[k]]
      }
      
      regime_ind[,k] <- (dat[, paste0('a',k)] == regime[,k])*1
      
      
    }
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')  
    
    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
  }
  return(regime_list)
}

# regimes <- list(list(0,c(0,0)),                
#                 list(0,c(0,1)),                
#                 list(0,c(1,0)),                
#                 list(0,c(1,1)),                
#                 list(1,c(0,0)),                
#                 list(1,c(0,1)),                
#                 list(1,c(1,0)),                
#                 list(1,c(1,1)))


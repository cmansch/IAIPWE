#### Make data function ####
# still 47.5
design2 <- function(vp, n, seed=1, s2=100, a1p=0.5, a2p=0.5, r2p = 0.5){
  set.seed(seed)
  # get the parameter values associated with the value pattern given by user
  if (vp==1){
    # to vary x higher, use 0.02, 3, 3
    vpbeta <- c(10, 0, 0, 0, 0.5, 12.5, 12.5, 0, 0)  # all equal
  } else if (vp==2) {
    vpbeta <- c(10, 0, 0, 5, 0.5, 12.5, 12.5, 0, 0)  # power one regime higher
  } else if (vp==3) {
    vpbeta <- c(15, -5, 0, 0, 0.5, 12.5, 12.5, 0, 0)  # make receiving A1=0 better 
  }
  # test vps
  # a1 <- c(0,0,0,0,1,1,1,1)
  # r2 <- c(0,0,1,1,0,0,1,1)
  # a2 <- c(0,1,0,1,0,1,0,1)
  # x11 <- 50
  # x12 <- 0.5
  # x21 <- 0.5
  # vs <- cbind(1, a1, a2, a1*a2, x11, x12, x21, a1*x11, a1*x12) %*% c(15, -5, 0, 0, 0.5, 12.5, 12.5, 0, 0)
  # mean(vs[c(1,3)])
  # mean(vs[c(1,4)])
  # mean(vs[c(2,3)])
  # mean(vs[c(2,4)])
  # mean(vs[c(5,7)])
  # mean(vs[c(5,8)])
  # mean(vs[c(6,7)])
  # mean(vs[c(6,8)])
  
  # randomly assign treatments to each individual
  a1 <- rbinom(n=n, size=1, prob=a1p)
  # create the error for individuals
  ei <- rnorm(n=n, mean=0, sd=sqrt(s2))
  # create the observed outcome based
  # aDesignMatrix <- cbind(1, a1, a2, a1*a2)
  
  # now the design specific for covariates (mean 0)
  x11 <- round(runif(n=n, min=25, max=75)) # mean 50 # age
  x12 <- rbinom(n=n, size=1, p=0.5) # mean 0.5 # gender
  x21 <- runif(n=n, min=0, max=1) # mean 0.5 # adherence
  # interim evaluated yi, using vpbeta[2] for the coefficient of a1
  # doing so correlates x21 with a1 and covariates
  # we comment this out, but use for the second simulation
  # x21 <- 10 + 0.02*x11 - 3*x12 + x13 + vpbeta[2]*a1
  # for covariate to be mean zero, we use c(0.02, -3, 1), do this to show
  # equivalence with previous approaches
  
  # r2 <- rbinom(n=n, size=1, prob=r2p)
  r2p <- exp(-1.1 + 0.02*x11 + .2*x12) / (1+exp(-1.1 + 0.02*x11 + .2*x12))
  r2 <- rbinom(n=n, size=1, prob=r2p)
  a2 <- rbinom(n=n, size=1, prob=a2p)  # this should be 0 or 1 for non-responders, and 0 for responders
  
  # y <- (aDesignMatrix %*% vpbeta) + cbind(x11, x12, x21) %*% c(0.02, 3, 3) + a1*x11 + a1*x12 + ei
  y <- (cbind(1, a1, a2, a1*a2, x11, x12, x21, a1*x11, a1*x12) %*% vpbeta) + ei
  
  # create the times of arrival for patients
  # these are done using the fixed times given to control for variance due
  # to enrollment process.
  # We multiple the process by 1000 to have the unit be 'days'.
  t1 <- round(runif(n, 0, 1), 3)*1000
  t2 <- t1 + 0.1*1000 + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  t3 <- t2 + 0.1*1000 + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  dat <- as.data.frame(cbind(a1, a2, y, r2, x11, x12, x21, t1, t2, t3))
  names(dat)[3] <- "y"
  # 
  # #### control arm data ####
  # nc <- floor(n/2)
  # 
  # x11c <- round(runif(n=nc, min=25, max=75)) # mean 50 # age
  # x12c <- rbinom(n=nc, size=1, p=0.5) # mean 0.5 # gender
  # eic <-  rnorm(n=nc, mean=0, sd=sqrt(s2))
  # 
  # yc <- 10 + 0.5*x11c + 25*x12c + eic    # need this to be 47.5
  # 
  # t1c <- round(runif(nc, 0, 1), 3)*1000
  # t2c <- t1c + 0.1*1000 + round(runif(nc, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # t3c <- t2c + 0.1*1000 + round(runif(nc, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  # datc <- as.data.frame(cbind(yc, x11c, x12c, t1c, t2c, t3c))
  # names(datc)[1] <- "y"
  
  return(dat)
}


#### Regime Indicator function ####
regimelist2 <- function(embRegimes, dat){
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


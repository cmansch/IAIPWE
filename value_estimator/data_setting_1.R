
design1 <- function(vp, n, seed=1, s2=100, a1p=0.5, a2p=0.5, r2p = 0.4){
  set.seed(seed)
  # get the parameter values associated with the value pattern given by user
  # cbind(1, a1, a2, a1*a2, x11, x12, x21, a1*x11, a1*x12)
  if (vp==1){
    # to vary x higher, use 0.02, 3, 3
    vpbeta <- c(10, 0, 0, 0, 0.5, 12.5, 12.5, 0, 0)  # all equal
  } else if (vp==2) {
    vpbeta <- c(10, 0, 0, 5, 0.5, 12.5, 12.5, 0, 0)  # power one regime higher
  } else if (vp==3) {
    vpbeta <- c(12.5, -2.5, 0, 0, 0.5, 12.5, 12.5, 0, 0)  # receiving A1=0 is better
  }
  
  ## to see what happens (0.4 is r2p and 0.5 a1p and a2p)
  # embr1 <- c(1, 0, 0, 0, 50, 0.5, 0.5, 0, 0)
  # embr2 <- c(1, 0, 1-0.4, 0, 50, 0.5, 0.5, 0, 0)
  # embr3 <- c(1, 1, 0, 0, 50, 0.5, 0.5, 0.5, 0.5)
  # embr4 <- c(1, 1, 1-0.4, 1-0.4, 50, 0.5, 0.5, 0.5, 0.5)
  # 
  # xm <- rbind(embr1, embr2, embr3, embr4)
  # 
  # xm %*% c(10, 0, 0, 0, 0.02, 3, 3, 0, 0)  # all equal
  # xm %*% c(10, 0, 0, 5, 0.02, 3, 3, 0, 0)  # power one regime higher
  # xm %*% c(10, 3, 5, 0, 0.02, 3, 3, 3, 3)  # make something else
  
  
  # randomly assign treatments to each individual
  a1 <- rbinom(n=n, size=1, prob=a1p)
  r2 <- rbinom(n=n, size=1, prob=r2p)
  a2 <- rbinom(n=n, size=1, prob=a2p)  # this should be 0 or 1 for non-responders, and 0 for responders
  a2 <- ifelse(r2 == 1, 0, a2)  # this ensures responders get 0, else 1. 
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
  
  return(dat)
}

regimelist <- function(embRegimes, dat, resptrt=list("r2" = list(0, 0, 0, 0))){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]])
  
  r <- 1
  for (regimen in embRegimes){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)
    
    for (k in 1:K){
      regime[,k] <- regimen[k]
      # update to reflect if a response at stage K is given
      if ((paste0("r", k) %in% colnames(dat)) & (!is.null(resptrt))) {
        regime[,k] <- regime[,k]*(1 - dat[, paste0('r',k)]) + 
          resptrt[[paste0('r',k)]][[r]] * dat[, paste0('r',k)]
      }
      regime_ind[,k] <- (dat[, paste0('a',k)] == regime[,k])*1
    }
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')
    
    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
    r <- r+1
  }
  
  return(regime_list)
}




# find chi-squared boundaries 
findChiSquareBound <- function(corr, alpha, mu=NULL, B=10000, tol=1e-4, 
                               startvalue=1, lambda=0.1, seed=1, K=2,
                               iota=c(1,1)){
  # iota=1 for pocock
  # iota=information proportion for OF
  # if equal information, then sqrt(1/0.5) = sqrt(2)
  set.seed(seed)
  # L regimes should be dim(corr)[1] / K
  L <- dim(corr)[1] / K

  Cmatrix <- cbind(diag(L-1), -1) 
  dfchi <- sum(svd(Cmatrix %*% corr[1:L,1:L] %*% t(Cmatrix))$d > 1e-8)

  if (is.null(mu)){
    mu <- rep(0, dim(corr)[1])   
  }
  
  Zb <- MASS::mvrnorm(n=B, mu=mu, Sigma=corr) # each row is a set of outcomes. 
  Zb <- matrix(Zb, nrow=B, byrow=FALSE)
  
  CmatrixBdiag <- Cmatrix
  if (K > 1){
    for (i in 2:K) {
      CmatrixBdiag <- bdiag(CmatrixBdiag, Cmatrix)
    }
  }
  
  CZ <- CmatrixBdiag %*% t(Zb) # these are the contrasts
  bigCorr <- CmatrixBdiag %*% corr %*% t(CmatrixBdiag)
  bigCorrInv <- MASS::ginv(as.matrix(bigCorr))
  
  chiscores <- matrix(0, nrow=B, ncol=K)
  for (i in 1:K){
    chiscores[,i] <- diag(t(Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),]) %*% 
                            MASS::ginv(as.matrix(bigCorr[((i-1) * (L-1)+1):(i*(L-1)), ((i-1) * (L-1)+1):(i*(L-1))])) %*% 
                            Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),])
  }

  calpha <- startvalue
  boundaries <- chiscores*0
  for(i in 1:K){
    boundaries[,i] <- calpha*iota[i]
  }
  
  
  # get the rejections
  # mean(rowSums(chiscores > boundaries) != 0)
  # as we increase startvalue, then the rejections should shrink
  # so if the starting rejections is below alpha, the start value is too high

  
  # alpha needs to be between 0 and 1
  stopifnot(alpha<1, alpha>0)
  # start with confidence level of zero, conf level is 1-alpha
  confLevel <- 0
  # initialize the bound
  bound <- startvalue
  
  # check if start value is too high
  confLevel <- 1 - mean(rowSums(chiscores > boundaries) != 0)
  
  
  if ((1-confLevel) < alpha){
    print('Starting value too high')
    while((1-confLevel) < alpha){
      calpha <- calpha - lambda
      
      for(i in 1:K){
        boundaries[,i] <- calpha*iota[i]
      }
      
      confLevel <- 1 - mean(rowSums(chiscores > boundaries) != 0)
    }
    print(paste0("New starting value: ", round(calpha,3)))
  }
  
  
  # begin loop to find true value. Requires the starting bound to achieve a
  # higher type I error rate than desired
  while(TRUE){
    while( TRUE){
      confLevel <- 1 - mean(rowSums(chiscores > boundaries) != 0)
      if ((1-confLevel) <= alpha){
        # if true, then the type 1 error is below alpha
        break
      } else{
        calpha <- calpha + lambda
        
        for(i in 1:K){
          boundaries[,i] <- calpha*iota[i]
        }
        # take step to higher calpha
      }
    }
    
    # check tolerance
    if (abs((1-confLevel ) - alpha) <= tol){
      # if true, we have reached a tolerance level we are happy with
      break
    }
    else{
      # update c to go 1*lambda prior to confLevel satisfied
      calpha <- calpha-lambda
      
      for(i in 1:K){
        boundaries[,i] <- calpha*iota[i]
      }
      # update lambda to be smaller
      lambda <- lambda/10
      # get confLevel at the new c
      confLevel <- 1 - mean(rowSums(chiscores > boundaries) != 0)
    }
  }
  return(list('bound' = calpha,
              'fnstatus' = confLevel))
}

# # Examples
# findChiSquareBound(corr=corr, alpha=0.05)
# findChiSquareBound(corr=corr, alpha=0.05, K=1, iota=1)

# # to check what the alpha is of a single test, use this
# pchisq(calpha, dfchi, ncp=0, lower.tail=FALSE, log.p=FALSE)
# pchisq(calpha*iota, dfchi, ncp=0, lower.tail=FALSE, log.p=FALSE)

findChiSquareSample <- function(corr, pow, cov, vp, calpha, interim=0.5,
                                lambda=10, N0=50, iota=c(1,1), seed=1, 
                                K=2, B=10000){
  stopifnot(interim<1, interim>0)
  stopifnot(lambda >= 1)
  stopifnot(max(abs(vp)) > 0)
  
  n2 <- N0
  n1 <- ceiling(n2*interim)
  
  # this is what changes as the sample changes
  # we use it to determine the alternative mean 
  mualt <- as.vector(vp * sqrt(c(rep(n1, dim(corr)[1]/2),rep(n2, dim(corr)[1]/2)) / cov))
  
  getrej <- function(mu){
    ## now we want to get the chi squares with the new mu
    set.seed(1)
    L <- dim(corr)[1] / K
    
    Cmatrix <- cbind(diag(L-1), -1) 
    dfchi <- sum(svd(Cmatrix %*% corr[1:L,1:L] %*% t(Cmatrix))$d > 1e-8)
    
    Zb <- MASS::mvrnorm(n=B, mu=mu, Sigma=corr) # each row is a set of outcomes. 
    Zb <- matrix(Zb, nrow=B, byrow=FALSE)
    
    CmatrixBdiag <- Cmatrix
    if (K > 1){
      for (i in 2:K) {
        CmatrixBdiag <- bdiag(CmatrixBdiag, Cmatrix)
      }
    }
    
    CZ <- CmatrixBdiag %*% t(Zb) # these are the contrasts
    bigCorr <- CmatrixBdiag %*% corr %*% t(CmatrixBdiag)
    bigCorrInv <- MASS::ginv(as.matrix(bigCorr))
    
    chiscores <- matrix(0, nrow=B, ncol=K)
    for (i in 1:K){
      chiscores[,i] <- diag(t(Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),]) %*% 
                              MASS::ginv(as.matrix(bigCorr[((i-1) * (L-1)+1):(i*(L-1)), ((i-1) * (L-1)+1):(i*(L-1))])) %*% 
                              Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),])
    }
    
    boundaries <- chiscores*0
    for(i in 1:K){
      boundaries[,i] <- calpha*iota[i]
    }
    # the expected number of rejections is
    rejH0 <- mean(rowSums(chiscores > boundaries) != 0)
    return(rejH0)
  }
  
  rejH0 <- getrej(mualt)
  
  if (rejH0 > pow){
    print('At the initial sample size, the desired power is attained')
    return(list('N' = c(n1, n2),
                'Power' = rejH0))
  }
  
  while(TRUE){
    while (TRUE){
      # in the case that the power is not enough yet, we want to increase the
      # sample size by lambda. Find new mean under alternative
      n2 <- n2 + lambda
      n1 <- ceiling(n2*interim)
      mualt <- as.vector(vp * sqrt(c(rep(n1, dim(corr)[1]/2),rep(n2, dim(corr)[1]/2)) / cov))
      # now we find power again
      rejH0 <- getrej(mualt)
      # check to see if we have achieved the power yet
      if (rejH0 > pow){
        break  # exit the loop
      }
    }
    # we reach this point when the power has been obtained.
    # now we want to go between the n2 -lambda and n2 by increments of 1
    n2 <- n2 - lambda
    n1 <- ceiling(n2*interim)
    while (TRUE){
      mualt <- as.vector(vp * sqrt(c(rep(n1, dim(corr)[1]/2),rep(n2, dim(corr)[1]/2)) / cov))
      # now we find power again
      rejH0 <- getrej(mualt)
      # once attained exit
      if (rejH0 > pow){
        return(list('N' = c(n1, n2),
                    'Power' = rejH0))
      } else{
        # otherwise increase n2 by 1
        n2 <- n2 + 1
        n1 <- ceiling(n2*interim)
      }
    }
  }
  
}

findChiSquareSampleFixed <- function(corr, pow, cov, vp, calpha,
                                lambda=10, N0=50, seed=1, 
                                K=1, B=10000){
  stopifnot(lambda >= 1)
  stopifnot(max(abs(vp)) > 0)
  
  n <- N0
  
  # this is what changes as the sample changes
  # we use it to determine the alternative mean 
  mualt <- as.vector(vp * sqrt(n / cov))
  
  getrej <- function(mu){
    ## now we want to get the chi squares with the new mu
    set.seed(1)
    L <- dim(corr)[1] / K
    
    Cmatrix <- cbind(diag(L-1), -1) 
    dfchi <- sum(svd(Cmatrix %*% corr[1:L,1:L] %*% t(Cmatrix))$d > 1e-8)
    
    Zb <- MASS::mvrnorm(n=B, mu=mu, Sigma=corr) # each row is a set of outcomes. 
    Zb <- matrix(Zb, nrow=B, byrow=FALSE)
    
    CmatrixBdiag <- Cmatrix
    if (K > 1){
      for (i in 2:K) {
        CmatrixBdiag <- bdiag(CmatrixBdiag, Cmatrix)
      }
    }
    
    CZ <- CmatrixBdiag %*% t(Zb) # these are the contrasts
    bigCorr <- CmatrixBdiag %*% corr %*% t(CmatrixBdiag)
    bigCorrInv <- MASS::ginv(as.matrix(bigCorr))
    
    chiscores <- matrix(0, nrow=B, ncol=K)
    for (i in 1:K){
      chiscores[,i] <- diag(t(Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),]) %*% 
                              MASS::ginv(as.matrix(bigCorr[((i-1) * (L-1)+1):(i*(L-1)), ((i-1) * (L-1)+1):(i*(L-1))])) %*% 
                              Cmatrix %*% t(Zb)[((i-1) * L + 1):(i*L),])
    }
    
    boundaries <- chiscores*0
    for(i in 1:K){
      boundaries[,i] <- calpha
    }
    # the expected number of rejections is
    rejH0 <- mean(rowSums(chiscores > boundaries) != 0)
    return(rejH0)
  }
  
  rejH0 <- getrej(mualt)
  
  if (rejH0 > pow){
    print('At the initial sample size, the desired power is attained')
    return(list('N' = n,
                'Power' = rejH0))
  }
  
  while(TRUE){
    while (TRUE){
      # in the case that the power is not enough yet, we want to increase the
      # sample size by lambda. Find new mean under alternative
      n <- n + lambda
      
      mualt <- as.vector(vp * sqrt(n / cov))
      # now we find power again
      rejH0 <- getrej(mualt)
      # check to see if we have achieved the power yet
      if (rejH0 > pow){
        break  # exit the loop
      }
    }
    # we reach this point when the power has been obtained.
    # now we want to go between the n2 -lambda and n2 by increments of 1
    n <- n - lambda
    
    while (TRUE){
      mualt <- as.vector(vp * sqrt(n / cov))
      # now we find power again
      rejH0 <- getrej(mualt)
      # once attained exit
      if (rejH0 > pow){
        return(list('N' = n,
                    'Power' = rejH0))
      } else{
        # otherwise increase n2 by 1
        n <- n + 1
      }
    }
  }
  
}


# bounds <- findChiSquareBound(corr=corr, alpha=0.05)
# findChiSquareSample(corr, 0.8, rep(300,8), c(4, 0, 0, 0, 0, 0, 0, 0), bounds$bound, interim=0.5,
#                                 lambda=10, N0=50, iota=c(1,1), seed=1, 
#                                 K=2, B=10000)
# bounds <- findChiSquareBound(corr=corr, alpha=0.05, K=1)
# findChiSquareSampleFixed(corr, pow, cov, vp, calpha,
#                                      lambda=10, N0=50, seed=1, 
#                                      K=1, B=10000)
#   
  
  
  
  
  
  
  



















# lambda is step size
# start value is the initialization
# tolerance is the tolerance window for 1-alpha
# alpha is the type I error rate
# corr is the correlation structure

findPocockBoundary <- function(corr, alpha, tol=0.000001, startvalue=0, lambda=0.01){
  # alpha needs to be between 0 and 1
  stopifnot(alpha<1, alpha>0)
  # start with confidence level of zero
  confLevel <- 0
  # initialize the bound
  bound <- startvalue

  # check if start value is too high
  confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                mean=rep(0, dim(corr)[1]), corr=corr)
  if ((1-confLevel) < alpha){
    print('Starting value too high')
    while((1-confLevel) < alpha){
      bound <- bound - lambda
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                    mean=rep(0, dim(corr)[1]), corr=corr)
    }
    print(paste0("New starting value: ", round(bound,3)))
  }



  # begin loop to find true value. Requires the starting bound to achieve a
  # higher type I error rate than desired
  while(TRUE){
    while( TRUE){
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                    mean=rep(0, dim(corr)[1]), corr=corr)
      if ((1-confLevel) <= alpha){
        # if true, then the type 1 error is below alpha
        break
      } else{
        bound <- bound+lambda  # take step to higher c
      }
    }

    # check tolerance
    if (abs((1-confLevel ) - alpha) <= tol){
      # if true, we have reached a tolerance level we are happy with
      break
    }
    else{
      # update c to go 1*lambda prior to confLevel satisfied
      bound <- bound-lambda
      # update lambda to be smaller
      lambda <- lambda/10
      # get confLevel at the new c
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                    mean=rep(0, dim(corr)[1]), corr=corr)
    }
  }
  return(list('bound' = bound,
              'fnstatus' = confLevel))
}

findOFBoundary <- function(corr, alpha, tol=0.000001, startvalue=0, lambda=0.01){
  # alpha needs to be between 0 and 1
  stopifnot(alpha<1, alpha>0)
  # start with confidence level of zero
  confLevel <- 0
  # initialize the bound
  bound <- startvalue

  # check if start value is too high
  upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
  confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                mean=rep(0, dim(corr)[1]), corr=corr)
  if ((1-confLevel) < alpha){
    print('Starting value too high')
    while((1-confLevel) < alpha){
      bound <- bound - lambda
      upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                    mean=rep(0, dim(corr)[1]), corr=corr)
    }
    print(paste0("New starting value: ", round(bound,3)))
  }



  # begin loop to find true value. Requires the starting bound to achieve a
  # higher type I error rate than desired
  while(TRUE){
    while( TRUE){
      upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                    mean=rep(0, dim(corr)[1]), corr=corr)
      if ((1-confLevel) <= alpha){
        # if true, then the type 1 error is below alpha
        break
      } else{
        bound <- bound+lambda  # take step to higher c
      }
    }

    # check tolerance
    if (abs((1-confLevel ) - alpha) <= tol){
      # if true, we have reached a tolerance level we are happy with
      break
    }
    else{
      # update c to go 1*lambda prior to confLevel satisfied
      bound <- bound-lambda
      # update lambda to be smaller
      lambda <- lambda/10
      # get confLevel at the new c
      upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
      confLevel <- mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                    mean=rep(0, dim(corr)[1]), corr=corr)
    }
  }
  return(list('bound' = bound,
              'fnstatus' = confLevel))
}



# mu under alternative hypothesis is the change to get sample size for power
# (value pattern repeated 2 times) * sqrt(c(n(s), n) ) / (sqrt(diag(V)))
# note that the variance will change depending on n(s)
#
# cor under the alternative should not change
findPocockSample <- function(corr, pow, cov, vp, bound, interim=0.5,
                             lambda=10, N0=50){
  # alpha needs to be between 0 and 1
  stopifnot(interim<1, interim>0)
  stopifnot(lambda >= 1)
  stopifnot(max(abs(vp)) > 0)

  n2 <- N0
  n1 <- ceiling(n2*interim)
  mualt <- as.vector(vp * sqrt(c(rep(n1, dim(corr)[1]/2),rep(n2, dim(corr)[1]/2)) / cov))

  # find the probability of not rejecting the null under alternative
  # rejH0 is 1-beta (i.e. the power)
  rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                mean=mualt, corr=corr)

  # if the designed power is above what we want, then we need to decrease n
  # this is unlikely since we start at n2 as 50
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
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                  mean=mualt, corr=corr)
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
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                  mean=mualt, corr=corr)
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

# note: cov should be a vector of the covariances of the value estimators
findOFSample <- function(corr, pow, cov, vp, bound, interim=0.5,
                             lambda=10, N0=50){
  # alpha needs to be between 0 and 1
  stopifnot(interim<1, interim>0)
  stopifnot(lambda >= 1)
  stopifnot(max(abs(vp)) > 0)

  n2 <- N0
  n1 <- ceiling(n2*interim)
  mualt <- as.vector(vp * sqrt(c(rep(n1, dim(corr)[1]/2),rep(n2, dim(corr)[1]/2)) / cov))

  # find the probability of not rejecting the null under alternative
  # rejH0 is 1-beta (i.e. the power)
  upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
  rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                              mean=mualt, corr=corr)

  # if the designed power is above what we want, then we need to decrease n
  # this is unlikely since we start at n2 as 50
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
      upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                  mean=mualt, corr=corr)
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
      upp <- c(sqrt(2)*rep(bound, dim(corr)[1]/2), rep(bound, dim(corr)[1]/2))
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=upp,
                                  mean=mualt, corr=corr)
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

# value pattern for findFixed should be for just the last analysis
findFixedSample <- function(corr, pow, cov, vp, bound,
                            lambda=10, N0=50){
  # alpha needs to be between 0 and 1
  stopifnot(lambda >= 1)
  stopifnot(max(abs(vp)) > 0)

  n <- N0
  mualt <- as.vector(vp * sqrt(n / cov))

  # find the probability of not rejecting the null under alternative
  # rejH0 is 1-beta (i.e. the power)
  rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                              mean=mualt, corr=corr)

  # if the designed power is above what we want, then we need to decrease n
  # this is unlikely since we start at n2 as 50
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
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                  mean=mualt, corr=corr)
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
      rejH0 <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=rep(bound, dim(corr)[1]),
                                  mean=mualt, corr=corr)
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

#
# # vp1cor <- as.matrix(read.csv(paste0(getwd(), "/sim_output/vp1_Zcor.csv"), row.names=1))
# vp1cor <- diag(1, 8)
# vp1cov <- diag(600,8)
# findPocockBoundary(corr=vp1cor, alpha=0.05, startvalue = 2)
# findOFBoundary(corr=vp1cor, alpha=0.05, startvalue = 2)
#
#
# findPocockSample(vp1cor, 0.8, vp1cov, vp=c(0,0,0,2,0,0,0,2),
#                  bound=2.48978)
#
# findOFSample(vp1cor, 0.8, vp1cov, vp=c(0,0,0,2,0,0,0,2),
#                  bound=2.255872)
#
# findPocockSample(vp1cor, 0.8, vp1cov, vp=c(0,2,2,0,0,2,2,0),
#                  bound=2.48978)
#
# findOFSample(vp1cor, 0.8, vp1cov, vp=c(0,2,2,0,0,2,2,0),
#              bound=2.255872)
#

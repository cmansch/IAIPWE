#### Make data function ####
# 750
design_test <- function(vp, n, seed=1, s2=100, a1p=0.5, a2p=0.5, r2p = 0.5){
  set.seed(seed)
  
  s2 <- 100
  x11 <- rnorm(n=n, mean=100, sd=sqrt(s2)) # baseline outcome measure
  x12 <- round(runif(n=n, min=25, max=75)) # mean 50 # age
  a1 <- rbinom(n=n, size=1, p=0.5) # mean 0.5
  
  x21 <- rnorm(n=n, mean=1.25 * x11, sd = 8^2) # interim outcome measure
  expit <- function(x){ exp(x) / (1 + exp(x))}
  r2 <- expit(-3 + 0.01*x11 + 0.02*x12 + 0.008*x21)
  # plot(density(r2)) # target 50% outcomes 
  a2 <- rbinom(n=n, size=1, prob=0.5) # a2 still 50 for everyone
  
  ei <- rnorm(n=n, mean=0, sd=500) # use here to control r^2 to around 60-70%
  # make the true outcome without any treatment 
  
  y <- 2.5*x11 + 5*x12 + 2*x21 + ei
  
  
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
regimelist_test <- function(embRegimes, dat){
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


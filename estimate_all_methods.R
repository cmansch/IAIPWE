#' The goal of this function is to take in a data set, function that identifes
#' regimes, the regimes themselves, and the outputs
#' 
#' 
#' t_s number of when (time wise) analyses should be done
#' if 0<t_s<1, assume quantile; else days
#' 
#' 
#' #' example call:
#' 
# for (file in list.files("./value_estimator", pattern=".R")){
#   source(paste0("./value_estimator/", file))
#   source("./smartim/R/vipw.R")
# }
# 
# library(doParallel)
# library(foreach)

unpack_IAIPW <- function(res){
  L <- length(res$values)
  n <- res$nus$ns
  endi <- dim(res$covariance)[2]
  starti <- endi - L + 1
  cov_ <- res$covariance[starti:endi, starti:endi]
  se_ <- sqrt(diag(cov_) / n)
  
  return( list('value' = res$values,
               'variance' = cov_,
               'se' = se_,
               'n' = n))
}

unpack_IPW <- function(res){
  L <- length(res$Values)
  n <- res$ns
  endi <- dim(res$Covariance)[2]
  starti <- endi - L + 1
  cov_ <- res$Covariance[starti:endi, starti:endi]
  se_ <- sqrt(diag(cov_) / n)
  
  return( list('value' = res$Values,
               'variance' = cov_,
               'se' = se_,
               'n' = n))
}

estimate_control <- function(dfc, t_s){
  maxt <- sum(unlist(lapply(colnames(dfc), function(x) substr(x,1,1))) == "t")
  if (0 < t_s & t_s < 1){
    t_s <- quantile(dfc[,paste0("t", maxt)], t_s)
  }
  if (isTRUE(t_s)) {# handle final analyses (t_s==TRUE)
    t_s <- max(dfc[,paste0("t", maxt)])
  }
  
  dft <- dfc[dfc[,paste0("t", maxt)] <= t_s,]
  
  return ( list('value' = mean(dft$y),
                'variance' = var(dft$y),
                'se' = sqrt( var(dft$y) / nrow(dft) ),
                'n' = nrow(dft) )
  )
  
}
# if t_s is TRUE, then we assume we are analyzing the entire data set. 
# if t_s is between 0 and 1 it is assumed to choose the proportion of completers for interim
# if t_s is >1 then it is assumed to be the actual time that analyses occur
estimateAllMethods <- function(df, regime_fn, regimes, q_list, pi_list, t_s, feasibleSetsIndicator){
  # handle if t_s is given as a quantile
  maxt <- sum(unlist(lapply(colnames(df), function(x) substr(x,1,1))) == "t")
  if (0 < t_s & t_s < 1){
    t_s <- quantile(df[,paste0("t", maxt)], t_s)
  }
  if (isTRUE(t_s)) {# handle final analyses (t_s==TRUE)
    t_s <- max(df[,paste0("t", maxt)])
  }
  
  # get the regimes of the whole data for interim AIPW estimator
  L <- length(regimes)
  regime_all <- regime_fn(embRegimes = regimes,
                           dat = df)
  # create reduced data set for AIPW, IPW estimators
  df_mid <- df[df[,paste0("t", maxt)] <= t_s,]
  regime_mid <- regime_fn(embRegimes = regimes,
                          dat = df_mid)
  
  # now get the estimates
  aa_res <- IAIPW(df=df, pi_list=pi_list, q_list=q_list, regime_all=regime_all, 
               feasibleSetsIndicator=feasibleSetsIndicator, t_s=t_s) 

  cd_res <- IAIPW(df=df_mid, pi_list=pi_list, q_list=q_list, regime_all=regime_mid, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=t_s) 
  
  ipw_res <- IAIPW(df=df_mid, pi_list=pi_list, q_list=NULL, regime_all=regime_mid, 
                   feasibleSetsIndicator=feasibleSetsIndicator, t_s=t_s) 
  
  # unpack results for easier use later
  return(list('aipw_aa' = unpack_IAIPW(aa_res),
              'aipw_cd' = unpack_IAIPW(cd_res),
              'ipw' = unpack_IAIPW(ipw_res)))
  
}


# q2 <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + I(1-a1):a2 +
#                                 I(1-a1):x21 + I(1-a1):a2:x21 + # for a1==0 arms
#                                 a1 + a1:x11 + a1:x12 + a1:a2 + # for a2==0 arms
#                                 a1:x21 + a1:a2:x21 - 1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# 
# q2r <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + I(1-a1):x21 +
#                                  a1 + x11:a1 + x12:a1 + x21:a1 - 1,
#                                solver.method = 'lm',
#                                predict.method = 'predict.lm')
# q1 <- modelObj::buildModelObj(model = ~ I(1-a1) + I(1-a1):x11 + I(1-a1):x12 + 
#                                 a1 + x11:a1 + x12:a1 - 1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# 
# q_list <- list(q1, list("r0" = q2,
#                         "r1" = q2r) )
# 
# p1 <- modelObj::buildModelObj(model = ~ 1,
#                               solver.method = 'glm',
#                               solver.args = list(family='binomial'),
#                               predict.method = 'predict.glm',
#                               predict.args = list(type='response'))
# 
# p2 <- modelObj::buildModelObj(model = ~ I(a1==0):I(r2==0) -1, # + I(a1==0):I(r2==1) + I(a1==1):I(r2==0),
#                               solver.method = 'glm',
#                               solver.args = list(family='binomial'),
#                               predict.method = 'predict.glm',
#                               predict.args = list(type='response'))
# pi_list <- list(p1, p2)
# 
# b<-1
# df <- design1(vp=1, n=500, seed=b*10000, s2=100, a1p=0.5, a2p=0.5, r2p = 0.4)
# regimes <- list(c(0,0), c(0,1), c(1,0), c(1,1))
# 
# estimateAllMethods(df, regimelist, regimes, q_list, pi_list, 700, TRUE)

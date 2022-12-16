#' This function is designed to estimate the variance of the value estimates
#' via simulation. We use the average asymptotic variance and to do such
#' In one sense, we could get an estimate of the asymptotic variance for N
#' large enough and use a single data set. Instead, we will use a large 
#' enough N and take enough MC samples to average the asymptotic variances
#' 
#' makeData should be a function to make the data sets we are interested in
#' vp
#' n
#' seed
#' loops
#' t_s_list is actually a vector
#' regimes is the encoding of regimes, its length should be L
#' 
#' This script will call functions in estimate_all_methods.R
library(doParallel)
library(foreach)
estimateVariances <- function(makeData, vp, n, regime_fn, regimes, loops, q_list, pi_list, 
                              t_s_list, feasibleSetsIndicator, cores=20,
                              delta,
                              out_location=NULL){
  
  registerDoParallel(cores=cores)
  
  L <- length(regimes)  # used later
  
  print('Simulating data')
  print(Sys.time())
  
  # estimate_all_methods.R is where these functions are used
  # to get a final analysis, use t_s_list = c(TRUE)
  
  # create an output list to save
  
  out_res <- foreach(n=rep(n, loops), loopno=c(1:loops), vp=rep(vp,loops),
                          regime_fn=replicate(n=loops, expr=regime_fn, simplify = FALSE),
                          regimes=replicate(n=loops, expr=regimes, simplify = FALSE),
                          q_list=replicate(n=loops, expr=q_list, simplify = FALSE),
                          pi_list=replicate(n=loops, expr=pi_list, simplify = FALSE),
                          t_s_list=replicate(n=loops, expr=t_s_list, simplify = FALSE),
                          feasibleSetsIndicator=rep(feasibleSetsIndicator, loops),
                          .packages=c("dplyr")
  ) %dopar% {
    # define the function that will execute in parallel
    # data 1 settings are 
    # vp, n, seed=1, s2=100, a1p=0.5, a2p=0.5, r2p = 0.4
    # require at least vp, n, and seed; other args may be passed
    
    # source functions needed here
    for (file in list.files("./value_estimator", pattern=".R")){
      source(paste0("./value_estimator/", file), local=TRUE)  # the local here keeps loop fast
    }
    source("./smartim/R/vipw.R", local=TRUE)
    source("./simulation_code/estimate_all_methods.R", local=TRUE)
    
    df <- makeData(vp=vp, n=n, seed=loopno * n)
    
    ## code added in case of a control arm 
    if(length(unique(df$a1) ) > 2){
      dfc <- df[df$a1 == 2,]  # create a data frame for the control arm
      df <- df[df$a1 != 2,]  # create a data frame for the smart
      cont_arm <- TRUE
    } else{
      cont_arm <- FALSE
    }
    ## end code update for control arm
    
    res <- vector("list", length(t_s_list)*2)
    for (i in 1:length(t_s_list)){
      t_s <- t_s_list[[i]]  # it needs to be a list for TRUE to be kepts
      res[[i]] <- estimateAllMethods(df, regime_fn, regimes, q_list, pi_list, t_s, feasibleSetsIndicator)
    }

    ## add code to estimate control arm 
    if(cont_arm){
      for (i in 1:length(t_s_list)){
        t_s <- t_s_list[[i]]  # it needs to be a list for TRUE to be kepts
        res[[i + length(t_s_list)]] <- estimate_control(dfc, t_s)
      }
    } else{
      for (i in 1:length(t_s_list)){
        t_s <- t_s_list[[i]]  # it needs to be a list for TRUE to be kepts
        res[[i + length(t_s_list)]] <- list('value' = 0,
                         'variance' = 0,
                         'se' = 0,
                         'n' = 0 )
      }
    }
    
    ## end code to estimate control arm
    names(res) <- c(paste0("t", 1:length(t_s_list)), paste0("c_t", 1:length(t_s_list)))
    
    res[["seed"]] <- loopno * n
    
    return(res)
  }
  stopImplicitCluster()
  rlist::list.save(out_res, file = paste0(out_location, "/varEstimateResults.Rdata"))
  
  print('Estimating the variance')
  print(Sys.time())
  
  # estimate the covariance and z-statistics from the data
  
  # initialize zs c() and vars matrix
  z_ipw <- c()
  z_aipw_cd <- c()
  z_aipw_aa <- c()
  
  var_ipw <- c()
  var_aipw_cd <- c()
  var_aipw_aa <- c()
  
  for (i in 1:(length(out_res[[1]]) - 1)){ # loops over ti
    # ipw
    z_temp <- lapply(out_res, function(x) (x[[paste0("t",i)]][['ipw']]$value - 
                                             x[[paste0("c_t",i)]]$value - delta) / 
                       sqrt( x[[paste0("t",i)]][['ipw']]$se^2 + 
                               x[[paste0("c_t",i)]]$se^2
                       ) 
                     )
    z_ipw <- cbind(z_ipw, 
                   matrix(unlist(z_temp), byrow=TRUE, nrow=length(z_temp))
                   )
    var_ipw <- c(var_ipw,
      # diag( Reduce("+", lapply(out_res, function(x) x[[paste0("t",i)]][['ipw']]$variance + 
      #                            x[[paste0("c_t",i)]]$variance )) / loops ) 
      diag( Reduce("+", lapply(out_res, function(x) 
        x[[paste0("t",i)]][['ipw']]$variance / (x[[paste0("t",i)]][['ipw']]$n / 
                                                  (x[[paste0("t",i)]][['ipw']]$n + x[[paste0("c_t",i)]]$n) ) + 
          x[[paste0("c_t",i)]]$variance / ( ifelse(test=x[[paste0("c_t",i)]]$n == 0, 
                                                   yes=1, no=x[[paste0("c_t",i)]]$n) / 
                                              (x[[paste0("t",i)]][['ipw']]$n + x[[paste0("c_t",i)]]$n) ) ) ) / loops ) 
    )
      
    # aipw_cd
    z_temp <- lapply(out_res, function(x) (x[[paste0("t",i)]][['aipw_cd']]$value - 
                                             x[[paste0("c_t",i)]]$value - delta) / 
                       sqrt( x[[paste0("t",i)]][['aipw_cd']]$se^2 + 
                               x[[paste0("c_t",i)]]$se^2
                       ) ) 
    z_aipw_cd <- cbind(z_aipw_cd, 
                       matrix(unlist(z_temp), byrow=TRUE, nrow=length(z_temp))
    )
    var_aipw_cd <- c(var_aipw_cd,
                     # diag( Reduce("+", lapply(out_res, function(x) 
                     #   x[[paste0("t",i)]][['aipw_cd']]$variance + 
                     #     x[[paste0("c_t",i)]]$variance )) / loops ) 
                     diag( Reduce("+", lapply(out_res, function(x) 
                       x[[paste0("t",i)]][['aipw_cd']]$variance / (x[[paste0("t",i)]][['aipw_cd']]$n / 
                                                                 (x[[paste0("t",i)]][['aipw_cd']]$n + x[[paste0("c_t",i)]]$n) ) + 
                         x[[paste0("c_t",i)]]$variance / ( ifelse(test=x[[paste0("c_t",i)]]$n == 0, 
                                                                  yes=1, no=x[[paste0("c_t",i)]]$n) / 
                                                             (x[[paste0("t",i)]][['aipw_cd']]$n + x[[paste0("c_t",i)]]$n) ) ) ) / loops ) 
    )
    # aipw_aa
    z_temp <- lapply(out_res, function(x) (x[[paste0("t",i)]][['aipw_aa']]$value - 
                                             x[[paste0("c_t",i)]]$value - delta) / 
                       sqrt( x[[paste0("t",i)]][['aipw_aa']]$se^2 + 
                               x[[paste0("c_t",i)]]$se^2
                       ) ) 
    z_aipw_aa <- cbind(z_aipw_aa, 
                       matrix(unlist(z_temp), byrow=TRUE, nrow=length(z_temp))
    )
    var_aipw_aa <- c(var_aipw_aa,
                         # diag( Reduce("+", lapply(out_res, function(x) 
                         #   x[[paste0("t",i)]][['aipw_aa']]$variance + 
                         #     x[[paste0("c_t",i)]]$variance )) / loops ) 
                     diag( Reduce("+", lapply(out_res, function(x) 
                       x[[paste0("t",i)]][['aipw_aa']]$variance / (x[[paste0("t",i)]][['aipw_aa']]$n / 
                                                                     (x[[paste0("t",i)]][['aipw_aa']]$n + x[[paste0("c_t",i)]]$n) ) + 
                         x[[paste0("c_t",i)]]$variance / ( ifelse(test=x[[paste0("c_t",i)]]$n == 0, 
                                                                  yes=1, no=x[[paste0("c_t",i)]]$n) / 
                                                             (x[[paste0("t",i)]][['aipw_aa']]$n + x[[paste0("c_t",i)]]$n) ) ) ) / loops )
    )
  } # end loop over ti
  
  # save vars, cov(z), corr(z) # save both z because we want diag to be 1, this is a good check.   
  write.csv(as.data.frame(cov(z_ipw)), paste0(out_location, "/Zcovariance", "ipw", ".csv"))
  write.csv(as.data.frame(cor(z_ipw)), paste0(out_location, "/Zcorrelation", "ipw", ".csv"))
  write.csv(as.data.frame(var_ipw), paste0(out_location, "/estimatedvariance", "ipw",".csv"))
  
  write.csv(as.data.frame(cov(z_aipw_cd)), paste0(out_location, "/Zcovariance", "cd", ".csv"))
  write.csv(as.data.frame(cor(z_aipw_cd)), paste0(out_location, "/Zcorrelation", "cd", ".csv"))
  write.csv(as.data.frame(var_aipw_cd), paste0(out_location, "/estimatedvariance", "cd",".csv"))
  
  write.csv(as.data.frame(cov(z_aipw_aa)), paste0(out_location, "/Zcovariance", "aa", ".csv"))
  write.csv(as.data.frame(cor(z_aipw_aa)), paste0(out_location, "/Zcorrelation", "aa", ".csv"))
  write.csv(as.data.frame(var_aipw_aa), paste0(out_location, "/estimatedvariance", "aa",".csv"))
  
  return()
} # end estimateVariances function


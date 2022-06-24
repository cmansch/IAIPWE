IAIPW <- function(df, pi_list, q_list, regime_all, feasibleSetsIndicator, t_s){
  K <- length(pi_list)
  n <- nrow(df)
  
  # estimate kappa
  kappa <- getKappa(df, t_s, K)
  
  # if no one has finished do not run
  stopifnot(max(kappa) == K+1)
  
  # estimate nu
  df[,'kappa'] <- kappa
  nus <- getNu(df, K)
  # estimate pis
  piFitted <- piFits(df, pi_list)
  pis <- piFitted[['ps']]
  p_fits <- piFitted[['p_fits']]
  # estimate values
  vHats<- estimate_values(df, q_list, regime_all, feasibleSetsIndicator, 
                              pis, nus)
  values <- vHats[['value']]
  dfs <- vHats[['df']]
  q_all <- vHats[['q_all']]
  
  # estimate the variance; these are already divided by n-p or n, respectively
  Bn <- getBn(df, p_fits, nus, q_all, dfs, values, feasibleSetsIndicator)
  An <- getAn(df, pis, p_fits, nus, q_all, values, regime_all, dfs, feasibleSetsIndicator)
  
  # An is a sparse Matrix object, it may be singular. 
  # In the event it is singular, take the ginv
  
  AnInv <- MASS::ginv(as.matrix(An))
  
  Vn <- AnInv %*% Bn %*% t(AnInv)
  # note that \widehat{\btheta} \sim N(0, Vn / n). Here n is actually nus$n
  # so the variance of \btheta is Vn / n
  
  # get all estimated parameters
  if (is.null(q_list)){
    estimatedParameters <- c(unlist(lapply(p_fits, function(x) x@fitObj$coefficients)), 
                             unlist(nus$nu),
                             values
    )
  } else{
    estimatedParameters <- c(unlist(lapply(p_fits, function(x) x@fitObj$coefficients)), 
                             unlist(nus$nu),
                             getQcoefs(q_all),
                             values
    )  
  }
  
  # unlist(lapply(q_all, function(y) unlist(lapply(y$q_fits, function(x) x@fitObj$coefficients)) ))
  
  return(list('values' = values,
              'covariance' = Vn,
              'params' = estimatedParameters,
              'An' = An,
              'Bn' = Bn,
              'nus' = nus,
              "q_all" = q_all,
              "regime_all" = regime_all,
              "dfs" = dfs))
}

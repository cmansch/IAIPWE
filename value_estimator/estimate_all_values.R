#' This function wraps the regression and value terms functions 
#' to estimate the value of the all regimes in list
#' It returns the estimated values, a list of the matrices with columns
#' corresponding to augmentation1-augmentation2K, IPW portion
#' and the list of the fit q-functions for each of the estimates

estimate_values <- function(df, q_list, regime_all, feasibleSetsIndicator, 
                            pis, nus){
  # get value estimate
  values <- c()
  vTermDfs <- list()
  q_all <- list()
  L <- length(regime_all)
  for (ell in 1:L){
    qs <- getQfits(df, q_list, regime_all[[ell]]$regime, feasibleSetsIndicator=feasibleSetsIndicator)
    vTerms <- getVterms(df, regime_all[[ell]]$regime_ind, pis, qs, nus)
    value <- sum(vTerms) / sum(df$kappa > 0)
    vTermDfs[[ell]] <- vTerms
    values <- c(values, value)
    q_all[[ell]] <- qs
  }
  return(list('value' = values,
              'df' = vTermDfs,
              'q_all' = q_all))
}




#' requires the data frame to have kappa stages in it. 
#' nu will be returned as a list with indexing opposite the paper
#' this index has k match stage k. 
#' nu[[K+1]] then is the probability that an individual has their final 
#' outcome observed given they are in the trial
getNu <- function(df, K){
  ns <- sum((df['kappa'] > 0))
  nu <- list()
  for (k in ((K+1):1)){
    nu[[k]] <- sum((df['kappa']>=k)*((df['kappa'] > 0))) / ns
  }
  # proportion of those who have received their last treatment,
  # but have not finished
  nd <- sum(df['kappa']>=(K+1)) / sum(df['kappa'] >= K)
  return(list('nu' = nu,
              'ns' = ns,
              'nd' = nd))
}

# get the Kappa(t_s) for individuals
# this is done by finding the total number of stages 
# an individual has reached based on arrival times t1,...,tK
getKappa <- function(df, t_s, K){
  return(rowSums(df[paste0('t', 1:(K+1))] <= t_s ))
}


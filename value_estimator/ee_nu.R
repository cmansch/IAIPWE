

getPsiNu <- function(df, nus){
  # psi nu for all k
  psi.nu <- c()
  for (k in (1:length(nus$nu))){
    psi.nu <- cbind(psi.nu, (df[,'kappa'] >= k) - (df[,'kappa'] >= 1)*nus$nu[[k]] )
  }
  # now for nud
  # psi.nu <- cbind(psi.nu, (df[,'kappa'] == (K+1)) - (df[,'kappa'] >= K)*nuFits$nd)
  return(psi.nu)  
}

getdPsiNu <- function(df, nus){
  # begin with all k
  # this will now be block diagonal. 
  psi.nu <- c()
  for (k in (1:length(nus$nu))){
    psi.nu <- c(psi.nu, -sum(df[,'kappa'] >= 1) )
  }
  # now for nud
  # psi.nu <- cbind(psi.nu, (df[,'kappa'] >= K) )
  return(diag(psi.nu))  # the base function diag makes the vector a diagonal matrix. 
}

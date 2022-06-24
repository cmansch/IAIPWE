


getBn <- function(df, p_fits, nus, q_all, dfs, values, feasibleSetsIndicator){
  # each function returns psi_i as a row with the EE changing wrt parameters across the columns
  # therefore, the sum psi_i psi_i^T can be written instead as the matrix X^TX for X has rows 
  # xi = psi_i. We then divide the estimator by n-p following Boos Stefanski pg 320. 
  # if we divide by n, that is also asymptotically unbiased. 
  
  # get psi pi
  psi.pi <- getPsiPi(df, p_fits)
  
  # get psi nu
  psi.nu <- getPsiNu(df, nus)
  
  # get psi beta for all regimes 
  # the getPsiBeta should do this for all of the embedded regimes 
  if (length(q_all) != 0){
    psi.beta <- getPsiBeta(df, q_all, feasibleSetsIndicator)
  } else{
    psi.beta <- NULL
  }
  
  # get values
  psi.v <- getPsiV(dfs, values, (df$kappa > 0))
  
  # psi <- as.data.frame(cbind(psi.pi1, psi.pi211, psi.pi201, psi.nu, psi.nu2, psi.nud, psibetas, psiv))
  psi <- as.data.frame(cbind(psi.pi, psi.nu, psi.beta, psi.v))
  
  # psi has dimension (n x p)
  
  # change made 3/7/2022: This should be n(s) not n(s) - p
  # return ( t(as.matrix(psi)) %*% as.matrix(psi) / (nus$ns - dim(psi)[2] ) )
  return ( t(as.matrix(psi)) %*% as.matrix(psi) / (nus$ns) )
  
}

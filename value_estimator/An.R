#' This code uses the derivatives from the estimating equations codes to get
#' the estimator for An, -dPsi/dTheta if EE are a vector (px1). 
#' We order the EE to match Bn, which is Pi1-PiK, nu1-nuK+1, ell1 Beta1-BetaK,
#' ellL Beta1-BetaK, V1-Vell. 
#' Based on section 7 of Boos Stefanski for Robust Sandwich Matrix

getAn <- function(df, pis, p_fits, nus, q_all, values, regime_all, dfs, feasibleSetsIndicator){
  
  # for pi, nu, and beta, these are the -derivatives of the estimating equations 
  # of those parameters with respect to only those parameters. This is because
  # they are 0 for all other entries.
  # however, dpsi.v is a matrix with each row as the -derivative of the regime-
  # specific estimating equation with respect to all parameters. 
  # so when we combine these terms, we create a block diagonal matrix 
  # from the first three (pi, nu, beta), then bind the rows of v to the bottom. 
  # we must bind a matrix of zeros to the columns of the (pi, nu, beta) because
  # it does not have the columns for dpsi.v
  
  # get dpsi pi, one call for all k and ell; this is NOT -dpsi.pi/dpi
  dpsi.pi <- getdPsiPi(df, p_fits)
  
  # get dpsi nu, one call for all k and ell; this is NOT -dpsi.nu/dnu
  dpsi.nu <- getdPsiNu(df, nus)
  
  # get psi beta for all regimes 
  # the getPsiBeta should do this for all of the embedded regimes 
  # the regime_all is added because now we deal with different design matrices
  # this is NOT -dpsi.beta/dbeta
  if (length(q_all) != 0){
    dpsi.beta <- getdPsiBeta(df, q_all, regime_all, feasibleSetsIndicator)
    # get values ee derivatives
    dpsi.v <- getdPsiV(df, pis, nus, regime_all, p_fits, q_all, feasibleSetsIndicator)
  } else{
    dpsi.v <- getdPsiV_ipw(df, pis, nus, regime_all, p_fits, q_all, feasibleSetsIndicator)
  }
  
  # create block diagonal matrix 
  if (length(q_all) != 0){
    dpsi.pi.nu.beta <- Matrix::bdiag(dpsi.pi, dpsi.nu, dpsi.beta)
  } else{
    dpsi.pi.nu.beta <- Matrix::bdiag(dpsi.pi, dpsi.nu)
  }
  
  # divide matrix by the sample size nus$ns, 
  return (rbind(cbind(dpsi.pi.nu.beta,
              matrix(0, nrow=nrow(dpsi.pi.nu.beta), ncol=nrow(dpsi.v)) 
              ),
              dpsi.v) / nus$ns
  )
  
}

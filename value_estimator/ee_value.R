

getPsiV <- function(dfs, values, cKappa){
  # cKappa indicates if individual observed yet (df$kappa > 0)
  makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  return ( t( makeNice(lapply(X=1:length(dfs), FUN=function(i) rowSums(dfs[[i]]) - values[i]*cKappa)
  )
  )
  )
}

# the challenge here is really getting the derivative with respect to all 
# parameters used in estimation. We recall the order need be pi, nu, betas, v
# we store the ee equations for individual portions of dPsiV in another script

getdPsiV <- function(df, pis, nus, regime_all, p_fits, q_all, feasibleSetsIndicator){
  # note that we need to get a row vector for each ell with the 
  # row having entries in the same order as Bn
  # i.e. Pi1-PiK, nu1-nuK+1, ell1 Beta1-BetaK, ellL Beta1-BetaK, V1-Vell. 
  # for ell != current, the betas should be zero. 
  # since V1-Vell is done as a block matrix, we can include this at the very end. 
  
  for (ell in 1:length(regime_all)){
    # start with pi
    dPsiV.dPi.ell <- c()
    dPsiV.dNu.ell <- c()
    dPsiV.dBeta.ell <- c()
    for (k in 1:length(p_fits)){
      dPsiV.dPi.ell <- cbind(dPsiV.dPi.ell, 
                             matrix(
                             getdPsiVdPi(df=df, pis=pis, nus=nus, 
                                         regime_ind=regime_all[[ell]]$regime_ind, 
                                         p_fits=p_fits, k=k, qs=q_all[[ell]]) 
                             ,nrow=1)
      )
    
      dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                             matrix(
                             getdPsiVdNu(df=df, pis=pis, nus=nus, 
                                         regime_ind=regime_all[[ell]]$regime_ind,
                                         k=k, qs=q_all[[ell]]) 
                             ,nrow=1)
      )
      if (k==length(p_fits)){
        dPsiV.dNu.ell <- cbind(dPsiV.dNu.ell, 
                               matrix(
                                 getdPsiVdNu(df=df, pis=pis, nus=nus, 
                                             regime_ind=regime_all[[ell]]$regime_ind,
                                             k=k+1, qs=q_all[[ell]]) 
                                 ,nrow=1)
        )
      }
      dPsiV.dBeta.ell <- cbind(dPsiV.dBeta.ell,
                               matrix(
                               getdPsiVdBeta(df=df, pis=pis, nus=nus, 
                                             regime_ind=regime_all[[ell]]$regime_ind,
                                             regime = regime_all[[ell]]$regime,
                                             q_fits = q_all[[ell]]$q_fits,
                                             k=k, feasibleSetsIndicator=feasibleSetsIndicator) 
                               ,nrow=1)
      )
                    
    }# end for loop over k stages
    if (ell==1){
      # initialize dPsiV as just that of Pi, Nu, Beta for regime 1
      dPsiV <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell, dPsiV.dBeta.ell)  
      zeroforbeta <- ncol(dPsiV.dBeta.ell)
    } else{
      # add a columns of zeros equal to the dimension of new betas to dPsiV
      dPsiV <- cbind(dPsiV, 
                     matrix(0, nrow=nrow(dPsiV), ncol=ncol(dPsiV.dBeta.ell)))
      dPsiV.ell <- cbind(dPsiV.dPi.ell, dPsiV.dNu.ell, matrix(0, ncol=zeroforbeta), 
                         dPsiV.dBeta.ell)
      # now the columns line up to bind the rows together
      dPsiV <- rbind(dPsiV, dPsiV.ell)
      # update the count of betas from previous regimes 
      zeroforbeta <- zeroforbeta + ncol(dPsiV.dBeta.ell)
    }

  }# end for loop over ell regime number 
  # now we should have the matrix of dPsiV / dPi dNu dBeta 11 - dBetaLK. 
  # add in the diagonal matrix corresponding to V
  dPsiV.V <- getdPsiVdV(nregimes=length(regime_all), nus=nus)
  
  dPsiV <- cbind(dPsiV, dPsiV.V)
  
  # this should return an entire row. It is still -dPsiV/dTheta
  # It is a row rather than a block because the estimating equations of V
  # are functions of the other parameters
  return (dPsiV) 
  
  
}

getPsiPi <- function(df, p_fits){
  g <- function(xi, thetapi){
    # xi is 1x4, thetapi is a vector of length 4
    return (exp(xi %*% thetapi) / (1 + exp(xi %*% thetapi)))
  }
  gprime <- function(xi, thetapi){
    # note that xi is 1x4 and the expfn is 1x1, result should be 4x1
    return (t(xi) %*% (exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^2)))
  }
  # second derivative of g
  gprime2 <- function(xi, thetapi){
    # xi is 1x4
    return (t(xi) %*% ( exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^3) ) %*% xi)
  }
  # used to format the output from lapply correctly
  # using t(sapply(...)) does not work if there is only 1 model parameter
  makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  
  
  psi.pi <- c()
  for (k in c(1:length(p_fits))){
    
    k_ind <- (df[,'kappa']>=k)*1
    dfkm <- model.matrix(p_fits[[k]]@modelObj@model, data=df)
    
    psii <- makeNice(
      lapply(X=1L:nrow(df), FUN = function(i) drop(df[,paste0('a',k)][i] - 
                                                     g(dfkm[i,,drop = FALSE], 
                                                       p_fits[[k]]@fitObj$coefficients)) *
               k_ind[i] *
               gprime(dfkm[i,,drop = FALSE], p_fits[[k]]@fitObj$coefficients)
      )
    )
    
    
    psi.pi <- cbind(psi.pi, psii) 
  }
  return(psi.pi)
}


getdPsiPi <- function(df, p_fits){
  g <- function(xi, thetapi){
    # xi is 1x4, thetapi is a vector of length 4
    return (exp(xi %*% thetapi) / (1 + exp(xi %*% thetapi)))
  }
  gprime <- function(xi, thetapi){
    # note that xi is 1x4 and the expfn is 1x1, result should be 4x1
    return (t(xi) %*% (exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^2)))
  }
  # second derivative of g
  gprime2 <- function(xi, thetapi){
    # xi is 1x4
    return (t(xi) %*% ( exp(xi %*% thetapi) / ((1 + exp(xi %*% thetapi))^3) ) %*% xi)
  }
  # used to format the output from lapply correctly
  # using t(sapply(...)) does not work if there is only 1 model parameter
  # makeNice <- function(...) {matrix(unlist(...), nrow=length(...), byrow=TRUE)}
  
  # general formula is given as eq 7.29 pg 319 Boos Stefanski
  # \sum_{i=1}^N gprime gprime^t - {(Y_i - g)g''}
  for (k in c(1:length(p_fits))){
    
    k_ind <- (df[,'kappa']>=k)*1
    dfkm <- model.matrix(p_fits[[k]]@modelObj@model, data=df)
    
    dpsi.pi_k <- Reduce("+", 
                        lapply(X=1L:nrow(df), FUN = function(i) k_ind[i] * 
                                 gprime(dfkm[i,,drop = FALSE], 
                                        p_fits[[k]]@fitObj$coefficients) %*% 
                                 t(gprime(dfkm[i,,drop = FALSE], 
                                          p_fits[[k]]@fitObj$coefficients)) 
                               # remove because expectation is mean 0, and adds computational burden. 
                               # - (df[i,'kappa']>=k) * {
                               #   drop(df[,paste0('a',k)][i] - 
                               #          g(model.matrix(p_fits[[k]]@modelObj@model, data=df[i,]), 
                               #            p_fits[[k]]@fitObj$coefficients)) * 
                               #     gprime2(model.matrix(p_fits[[k]]@modelObj@model, data=df[i,]), p_fits[[k]]@fitObj$coefficients)
                               # }
                        )
    )
    # this should be block diagonal
    if (k == 1){
      dpsi.pi <- dpsi.pi_k
    } else{
      dpsi.pi <- Matrix::bdiag(dpsi.pi, dpsi.pi_k)
    }
  }
  return(-dpsi.pi)
  return(-dpsi.pi)
}


    


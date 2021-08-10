#' This file contains scripts used to autocalculate regimes under certain
#' scenarios. These functions are not intended for general use but may be used
#' in the case of the test functions and the data generation found in the
#' accompanying paper
#'
#' This function takes an embedded regime and creates a treatment list
#' as if the regime was followed. It also creates an indicator if the
#' embedded regime matches the regime observed.
#'
#' @param embRegimes list,
#' @param dat data frame, must include data matching the regimes, the names
#'            for treatment should be a1.
#' @export

regimelist <- function(embRegimes, dat){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]])

  r <- 1
  for (regimen in embRegimes){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)

    for (k in 1:K){
      regime[,k] <- regimen[k]
      regime_ind[,k] <- (dat[, paste0('a',k)] == regimen[k])*1
    }
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')

    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
    r <- r+1
  }

  return(regime_list)
}

# generate regimes
# response generator 

c(0,0,0)  # give a1,a2,responders get 0


regimeListResponse <- function(embRegimes, dat, nrespinds){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]]) - nrespinds
  for (regimen in embRegimes){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)
    for (k in 1:K){
      regime[,k] <- regimen[k]
      regime_ind[,k] <- (dat[, paste0('a',k)] == regimen[k])*1
      # check against the response status
    }
    for (k in (K-nrespinds):K){
      
    }
  }
  
  
}




regimeListResponse <- function(embRegimes, dat){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]])
  
  r <- 1
  for (regimen in embRegimes){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)
    
    for (k in 1:K){
      regime[,k] <- regimen[k]
      regime_ind[,k] <- (dat[, paste0('a',k)] == regimen[k])*1
      # check against the response status
    }
    
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')
    
    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
    r <- r+1
  }
  
  return(regime_list)
}

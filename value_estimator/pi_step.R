#' propensity model building code giving here
#' we assume a list of propensity models is passed to the function
#' this list of models can be tabular outcomes through use of 
#' binomial functions with the terms in the model being
#' indicating up to the node where randomization occurs
#' here the regime should be the columns of df corresponding to the actual 
#' treatment received

pstep <- function (pmodel, data, response, k) {
  # fit the qstep model
  pfit <- modelObj::fit(object = pmodel, data = data, response = response)
  # get unmodified predicted values for newdata
  pk <- modelObj::predict(object = pfit, newdata = data)
  # these are now the probability that the outcome is 1
  # however, we want the probability that individual receives the 
  # treatment that they did indeed get
  # ps[[k]] <- {pk*{regime[,k] == 1L} + {1.0 - pk}*{regime[,k] == 0L}}
  ps <- {pk*{ data[,paste0('a',k)] == 1L} + {1.0 - pk}*{data[,paste0('a',k)] == 0L}}
  return (list('pk' = pk, 
               'ps' = ps,
               'pfit' = pfit))
}

piFits <- function(df, p_list){  # , regime){
  K <- length(x = p_list)
  fits <- list()
  ps <- data.frame(matrix(ncol=K, nrow=nrow(df)))
  # we want the names of the columns to be pi1,...,piK. 
  colnames(ps) <- paste0("pi", c(1:K))
  
  for (k in 1:K){
    # only estimate for individuals with observed ak
    ones <- (df['kappa']>=k)  
    
    # fits[[k]] <- modelObj::fit(object = p_list[[k]], 
    #                            data = df, 
    #                            response = df[[paste0('a',k)]])
    # pk <- modelObj::predict(object = fits[[k]],
    #                         data = df)
    # # these are now the probability that the outcome is 1
    # # however, we want the probability that individual receives the 
    # # treatment that they did indeed get
    # # ps[[k]] <- {pk*{regime[,k] == 1L} + {1.0 - pk}*{regime[,k] == 0L}}
    # ps[[k]] <- {pk*{ df[,paste0('a',k)] == 1L} + {1.0 - pk}*{df[,paste0('a',k)] == 0L}}
    
    pks <- pstep(pmodel = p_list[[k]],
                data = df[ones,, drop = FALSE],
                response = df[ones, paste0('a',k), drop = FALSE],
                k = k)
    
    fits[[k]] <- pks$pfit
    ps[,k] <- 99  # give everyone something. 
    # these numbers do not factor in because the estimates pi should always be multiplied by the 
    # indicator for kappa to ensure that it is 0 for anyone who has not reached
    # stage k. This is factored into the estimating equations. 
    ps[ones, k] <- pks$ps  # only update those with estimates
  }

  return (list('ps' = ps,
               'p_fits' = fits))
  
}

# p1 <- modelObj::buildModelObj(model = ~ 1,
#                               solver.method = 'glm',
#                               solver.args = list(family='binomial'),
#                               predict.method = 'predict.glm',
#                               predict.args = list(type='response'))
# 
# p2 <- modelObj::buildModelObj(model = ~ I(a1==0):I(r2==0) + I(a1==0):I(r2==1) + I(a1==1)*I(r2==0),
#                               solver.method = 'glm',
#                               solver.args = list(family='binomial'),
#                               predict.method = 'predict.glm',
#                               predict.args = list(type='response'))
# 
# 
# a1 <- rbinom(n=100, size=1, prob=0.5)
# r2 <- rbinom(n=100, size=1, prob=0.5)
# a2 <- rbinom(n=100, size=1, prob=0.5)
# a2 <- ifelse(a1*r2==1, 1, a2)
# df <- data.frame(a1=a1, a2=a2, r2=r2)
# 
# pis <- piFits(df, p_list=list(p1, p2))

# this should only use individuals at stage k or higher. 


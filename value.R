#' This file is used to keep the sub routines for the AIPW estimator
#' To install package on unix, use the following commands
#' R CMD build smartim
#' R CMD INSTALL smartim
library(dplyr)

#### get kappa for time ti ####
getKappa <- function(df, ti, maxEvalNumber){
  return(rowSums(df[paste0('t', 1:maxEvalNumber)] <= ti ))
}
#### regressions ####
# get the model fits for a given regime
qstep <- function (qmodel, data, response, newdata, regime, txName) {
  # fit the qstep model
  qfit <- modelObj::fit(object = qmodel, data = data, response = response)
  # get unmodified predicted values for newdata
  hats_unmod <- modelObj::predict(object = qfit, newdata = newdata)
  # set tx to recommended
  newdata[,txName] <- regime
  # predict the newdata with regime consistent estimates
  hats_mod <- modelObj::predict(object = qfit, newdata = newdata)
  return (list('hats_mod' = hats_mod,
               'hats_unmod' = hats_unmod,
               'qfit' = qfit))
}
getQfits <- function(df, q_list, regime, feasibleSetsIndicator=NULL){
  # total number of decision points
  K <- length(q_list)
  # create storage for modified and unmodified predictors
  unmod_regime_vhats <- matrix(data = 0.0, nrow = nrow(df), ncol = K + 1)
  mod_regime_vhats <- matrix(data = 0.0, nrow = nrow(df), ncol = K + 1)

  # for Kappa = decisions + 1 keep the response
  unmod_regime_vhats[,K + 1] <- df$y
  mod_regime_vhats[,K + 1] <- df$y

  q_fits <- list()

  for (k in (K:1)){
    # fit the model using only those with multiple treatments
    if (is.null(feasibleSetsIndicator)){
      # only fit those who have finished the trial and were not responders
      # because responders do not have multiple treatment options
      ones <- (df['kappa']>=(k+1)) # those who have finished
      new <- (df['kappa']>=k) # those who can attain a prediction from model
      # get the fitted models and the predicted values
      vk <- qstep(qmodel = q_list[[k]],
                  data = df[ones,, drop = FALSE],
                  response = mod_regime_vhats[ones, k+1, drop = FALSE],
                  newdata = df[new,, drop = FALSE],
                  regime = regime[new,k],
                  txName = paste0('a',k) )
      # update
      mod_regime_vhats[, k] <- mod_regime_vhats[,k+1]
      mod_regime_vhats[new, k] <- vk$hats_mod
      unmod_regime_vhats[, k] <- unmod_regime_vhats[,k+1]
      unmod_regime_vhats[new, k] <- vk$hats_unmod
      q_fits[[k]] <- vk$qfit

    } else {
      #### This will need to be fixed
      print('There are multiple feasible sets')
      # ones <- (df['kappa']==3) & (df['r2'] == 0)  # make sure this is boolean
      # new2 <- (df['kappa']>=2) & (df['r2'] == 0)  # make sure this is boolean
    }
  }

  colnames(mod_regime_vhats) <- paste0('q', 1:(K+1))
  colnames(unmod_regime_vhats) <- paste0(paste0('q', 1:(K+1)), '_nochange')
  return(list('q_fits' = q_fits,
              'mod_regime_vhats' = mod_regime_vhats,
              'unmod_regime_vhats' = unmod_regime_vhats))

}

#### get the nu values ####
getNu <- function(df, K){
  ns <- sum((df['kappa'] > 0))
  nu <- list()
  for (k in ((K+1):1)){
    nu[[k]] <- sum((df['kappa']>=k)*((df['kappa'] > 0))) / ns
  }
  return(list('nu' = nu,
              'ns' = ns))
}

#### get pi ####
getPi <- function(df){
  # start at time step one as each consecutive table requires the
  # n from the previous step
  pi1_table <- subset(df, kappa >= 1) %>% count(a1) %>% mutate(freq=n/sum(n))
  df <- df %>% left_join(pi1_table[,c('a1', 'freq')], by="a1")
  names(df)[names(df) == 'freq'] <- 'pi1'
  names(pi1_table) <- c('a1', 'n1', 'freq1')

  pi2_table <- subset(df, kappa >= 2) %>%
    count(a1, a2) %>%
    group_by(a1) %>%
    mutate(freq=n/sum(n))
  names(pi2_table) <- c('a1', 'a2', 'n2', 'pi2')
  df <- df %>% left_join(pi2_table[,c('a1', 'a2', 'pi2')], by=c("a1", "a2"))

  # zero out pi1, pi2 if kappa=0
  df[df['kappa'] == 0,c('pi1', 'pi2')] <- 0

  # zero out pi2 if kappa=1
  df[df['kappa'] == 1,c('pi2')] <- 0

  return(df[,c('pi1','pi2')])
}

#### get lead, augmentation1, ..., augmentationk ####
getVterms <- function(df, steps, nu){

  # for (k in 1:(steps+1)){
  # start with kappa=1
  dfk1 <- df[df['kappa'] == 1,]
  if (nrow(dfk1) > 0){
    dfk1['lead'] <- 0
    dfk1['aug1'] <- (dfk1['a1_ind'] - dfk1['pi1'])*dfk1['q1'] /
      (dfk1['pi1'])
    dfk1['aug2'] <- 0
  }
  # now do kappa=2
  dfk2 <- df[df['kappa'] == 2,]
  if (nrow(dfk2) > 0){
    dfk2['lead'] <- 0
    dfk2['aug1'] <- (dfk2['a1_ind'] - dfk2['pi1'])*dfk2['q1'] /
      (dfk2['pi1'])
    dfk2['aug2'] <- dfk2['a1_ind']*(dfk2['a2_ind'] - dfk2['pi2'])*
      dfk2['q2'] /
      (nu[[2]]*dfk2['pi1']*dfk2['pi2'])
  }
  # now do kappa=3, this means someone has finished
  dfk3 <- df[df['kappa'] == 3,]
  if (nrow(dfk3) > 0){
    dfk3['lead'] <- (dfk3['a1_ind']*dfk3['a2_ind'])*dfk3['y'] / (nu[[3]]*dfk3['pi1']*dfk3['pi2'])
    dfk3['aug1'] <- (dfk3['a1_ind'] - dfk3['pi1'])*dfk3['q1'] /
      (dfk3['pi1'])
    dfk3['aug2'] <- dfk3['a1_ind']*(dfk3['a2_ind'] - dfk3['pi2'])*
      dfk3['q2'] /
      (nu[[2]]*dfk3['pi1']*dfk3['pi2'])
  }
  # now stack
  # stacking empty df does not affect this
  # this df_order will now be only of size ns
  df_order <- rbind(dfk1, dfk2, dfk3)
  # }
  return(df_order)
}

#### calculate An and Bn for all regimes ####
#### get psi beta for a particular regime

getPsiBeta <- function(df_order, q2, q1){
  q2_design <- model.matrix(q2@modelObj@model, data=df_order)
  q1_design <- model.matrix(q1@modelObj@model, data=df_order)
  # now we make the design matrix for the regime consistent estimates
  # these are done from the last decision forward since we only modify the action
  # of the same index as the regression fit
  df_mod <- df_order
  df_mod[,'a2'] <- df_order[,'regime2']
  q2_design_mod <- model.matrix(q2@modelObj@model, data=df_mod)
  df_mod[,'a1'] <- df_order[,'regime1']
  q1_design_mod <- model.matrix(q1@modelObj@model, data=df_mod)
  psi.beta1 <- diag((df_order[,'q2'] - df_order[,'q1_nochange'])*(df_order[,'kappa'] >= 2)) %*%
    q1_design
  # beta2 # should be a matrix for each entry
  psi.beta2 <- diag((df_order[,'y'] - df_order[,'q2_nochange'])*(df_order[,'kappa'] == 3)) %*%
    q2_design
  return(cbind(psi.beta1, psi.beta2))
}


#### Bn ####
getBn <- function(dfList, qFitsList, nuFits, vList){
  # these contain the different results from each.
  # pi, nu shared, these values are shared across common df
  df_order <- dfList[[1]]  # just use the first one for pi and nu
  ## psi values
  # these may need to be separated out for each a2,h2 combination
  # in order to get pi1 the parameter for this, we use the df_order value.
  psi.pi1 <- (df_order[,'a1'] == 1)*(df_order[,'kappa'] >= 1) -
    (df_order[,'kappa'] >= 1)*abs(df_order[,'pi1'] - 1*(df_order[,'a1'] == 0))
  # change made 6.9 based on how this is actually estimated, needs a1 as either 0 or 1
  psi.pi211 <- (df_order[,'a1'] == 1)*((df_order[,'kappa'] >= 2)*(df_order[,'a2'] == 1) - (df_order[,'kappa'] >= 2)*
                  (abs(df_order[,'pi2'] - (df_order[,'a2'] == 0) )*(df_order[,'a1'] == 1)*(df_order[,'kappa'] >=2)))
  psi.pi201 <- (df_order[,'a1'] == 0)*((df_order[,'kappa'] >= 2)*(df_order[,'a2'] == 1) - (df_order[,'kappa'] >= 2)*
                                         (abs(df_order[,'pi2'] - (df_order[,'a2'] == 0) )*(df_order[,'a1'] == 0)*(df_order[,'kappa'] >=2)))
  # nu values, at the final time point these will all be 0
  psi.nu <- (df_order[,'kappa'] == 3) - (df_order[,'kappa'] >= 1)*nuFits$nu[[3]]
  # psi.nu1 <- (df_order[,'kappa'] >= 1)*(df_order[,'kappa'] >= 1) - (df_order[,'kappa'] >= 1)*nuFits$nu[[1]]
  psi.nu2 <- (df_order[,'kappa'] >= 2) - (df_order[,'kappa'] >= 1)*nuFits$nu[[2]]

  nRegimes <- length(vList)
  psibetas <- c()
  psiv <- c()
  for (nreg in 1:nRegimes){
    # now we must get psi.beta for each of the regimes considered
    psibetas <- cbind(psibetas, getPsiBeta(dfList[[nreg]], q2=qFitsList[[nreg]]$q_fits[[2]],
                                           q1=qFitsList[[nreg]]$q_fits[[1]]))
    # we also must get the value estimators for each regime considered
    vi <- (dfList[[nreg]]['lead'] - dfList[[nreg]]['aug1'] - dfList[[nreg]]['aug2']) / nuFits$ns
    psiv <- cbind(psiv, (vi*nuFits$ns - vList[[nreg]])[,1])
  }
  colnames(psiv) <- paste0('v',1:(length(vList)) )
  # psi <- as.data.frame(cbind(psi.pi1, psi.pi211, psi.pi201, psi.nu, psi.nu1, psi.nu2, psibetas, psiv))
  psi <- as.data.frame(cbind(psi.pi1, psi.pi211, psi.pi201, psi.nu, psi.nu2, psibetas, psiv))

  Bn <- (t(as.matrix(psi)) %*% as.matrix(psi)) / dim(psi)[1]
  return(Bn)
}

#### An ####
getBetaDims <- function(qFitsList){
  # the first layer is the embedded regime
  # the second layer is which step is being fit
  # the third layer is the dimension of the coefficients
  betaDimsOut <- list()
  singleRegime <- list()
  for (i in 1:length(qFitsList)){
    for (j in 1:length(qFitsList[[i]]$q_fits)){
      singleRegime[[j]] <- length(qFitsList[[i]]$q_fits[[j]]@fitObj$coefficients)
    }
    betaDimsOut[[i]] <- singleRegime
    singleRegime <- list()  # reset for next regime
  }
  return(betaDimsOut)
}

getdPsiBeta <- function(df_order, q2, q1, dimbeta1,
                        dimbeta2, leadingzeros, dimBn, i){
  q2_design <- model.matrix(q2@modelObj@model, data=df_order)
  q1_design <- model.matrix(q1@modelObj@model, data=df_order)
  # now we make the design matrix for the regime consistent estimates
  # these are done from the last decision forward since we only modify the action
  # of the same index as the regression fit
  df_mod <- df_order
  df_mod[,'a2'] <- df_order[,'regime2']
  q2_design_mod <- model.matrix(q2@modelObj@model, data=df_mod)
  df_mod[,'a1'] <- df_order[,'regime1']
  q1_design_mod <- model.matrix(q1@modelObj@model, data=df_mod)

  # matrix of zeroes dimbeta1 by leadingzeros
  # then dpsibeta1 dbeta1, dpsibeta1 dbeta2
  # then matrix of zeroes to pad fill through dimension
  dpsibeta1 <- cbind( matrix(0, nrow=dimbeta1, ncol=leadingzeros) ,
                      (df_order[,'kappa'] >= 2)[i] * q1_design[i,] %*%
                        t(q1_design[i,]),
                      -(df_order[,'kappa'] >= 2)[i] * q1_design[i,] %*%
                        t(q2_design_mod[i,]),
                      matrix(0, nrow=dimbeta1, ncol=dimBn - leadingzeros - dimbeta1 - dimbeta2)
  )
  #beta2
  dpsibeta2 <- cbind( matrix(0, nrow=dimbeta2, ncol=leadingzeros+dimbeta1) ,
                      (df_order[,'kappa'] == 3)[i] *
                        q2_design[i,] %*% t(q2_design[i,]),
                      matrix(0, nrow=dimbeta2, ncol=dimBn - leadingzeros - dimbeta1 - dimbeta2) )

  return(rbind(dpsibeta1, dpsibeta2))
}


repnull <- function(n, times){
  if(is.null(times)){
    return(NULL)
  } else{
    return(rep(n, times))
  }
}


getPartialValueSingleRegime <- function(df_order, q2, q1, nuFits,
                                        leadbetas, followbetas,
                                        leadv, followv, i){
  # note that lead and follow should be NULL for first and last
  # regimes respectively to avoid including an extra 0
  q2_design <- model.matrix(q2@modelObj@model, data=df_order)
  q1_design <- model.matrix(q1@modelObj@model, data=df_order)
  # now we make the design matrix for the regime consistent estimates
  # these are done from the last decision forward since we only modify the action
  # of the same index as the regression fit
  df_mod <- df_order
  df_mod[,'a2'] <- df_order[,'regime2']
  q2_design_mod <- model.matrix(q2@modelObj@model, data=df_mod)
  df_mod[,'a1'] <- df_order[,'regime1']
  q1_design_mod <- model.matrix(q1@modelObj@model, data=df_mod)

  nu <- nuFits$nu[[3]]
  nu2 <- nuFits$nu[[2]]

  dv <- -(df_order[i,'kappa'] >= 1)*
    c((-df_order[i,'lead']/df_order[i,'pi1'] +
        (df_order[i,'a1']==df_order[i,'regime1'])*df_order[i,'q1'] / (df_order[i,'pi1']^2) +
        df_order[i,'aug2']/df_order[i,'pi1'])*(2*(df_order[i,'a1'] == 1)-1)  # pi1
      ,(-df_order[i,'lead']/df_order[i,'pi2'] -
        -((df_order[i,'kappa'] >= 2)/nu2)*((df_order[i,'a1']==df_order[i,'regime1'])/df_order[i,'pi1'])*
        ((df_order[i,'a2']==df_order[i,'regime2'])*df_order[i,'q2'] / (df_order[i,'pi2']^2)))*
        (2*(df_order[i,'a1'] == 1)*(df_order[i,'a2'] == 1)-1)#pi21
      ,(-df_order[i,'lead']/df_order[i,'pi2'] -
        -((df_order[i,'kappa'] >= 2)/nu2)*((df_order[i,'a1']==df_order[i,'regime1'])/df_order[i,'pi1'])*
        ((df_order[i,'a2']==df_order[i,'regime2'])*df_order[i,'q2'] / (df_order[i,'pi2']^2)))*
        (2*(df_order[i,'a1'] == 0)*(df_order[i,'a2'] == 1)-1)#pi20
      ,-df_order[i,'lead']/nu  #nu
      # ,0  #nu1  this is not needed
      ,--df_order[i,'aug2']/nu2   #nu2
      ,repnull(0,leadbetas)
      ,(((df_order[i,'a1']==df_order[i,'regime1']) - (df_order[i,'pi1'])) / (df_order[i,'pi1']))*
        q1_design_mod[i,] #beta1
      ,((df_order[i,'kappa'] >= 2)/nu2) * ((df_order[i,'a1']==df_order[i,'regime1']) / (df_order[i,'pi1'])) *
        (((df_order[i,'a2']==df_order[i,'regime2']) - (df_order[i,'pi2'])) / (df_order[i,'pi2'])) *
        q2_design_mod[i,]  #beta2
      ,repnull(0,followbetas)
      ,repnull(0,leadv)  # start v
      , -1
      ,repnull(0,followv) )
  return(dv)

}


getAn <- function(dimBn, vList, K, dfValueList,
                  qFitsList, nuFits){
  An <- matrix(0, nrow=dimBn, ncol=dimBn)
  # create the derivative of psi wrt parameters in same order as Bn
  # so for each i, we will need to calculate the different parts
  # we first need to get the dimensions of the parameters
  dimv <- length(vList)
  dimpi <- K
  dimnu <- K+1
  # now get dim betas
  # we will order the list in a layered manner
  dimbetas <- getBetaDims(qFitsList)
  # now we get the derivatives of the estimating equations
  # by row (i.e. fix psi, get all derivatives going left to right)
  # then we will rbind the rows to get An
  df <- dfValueList[[1]]
  for (i in 1:nrow(df)){
    # dpsi1,2 (pi) dtheta
    dpsi1dtheta <- c((df[,'kappa'] >= 1)[i]*1, rep(0,dimBn-1))
    # change made on 6.9 to match the paper.
    dpsi2dtheta1 <- c(rep(0,1),(df[,'kappa'] >= 2)[i]*(df[i,'a1']==1)*1, rep(0,dimBn-2))
    dpsi2dtheta2 <- c(rep(0,2),(df[,'kappa'] >= 2)[i]*(df[i,'a1']==0)*1, rep(0,dimBn-3))
    AnTemp <- rbind(dpsi1dtheta, dpsi2dtheta1, dpsi2dtheta2)

    # dpsi3,4,5 (nu) dtheta
    dnu <- c(rep(0,3), (df[,'kappa'] >= 1)[i], rep(0, dimBn-4))
    # dnu1 <- c(rep(0,4), (df[,'kappa'] >= 1)[i], rep(0, dimBn-5))
    dnu2 <- c(rep(0,4), (df[,'kappa'] >= 1)[i], rep(0, dimBn-5))
    # AnTemp <- rbind(AnTemp, dnu, dnu1, dnu2)
    AnTemp <- rbind(AnTemp, dnu, dnu2)

    leadingzeros <- nrow(AnTemp)
    leadbetas <- 0
    leadv <- 0
    followbetas <- sum(unlist(dimbetas))  # total parameters for regression
    followv <- length(vList)

    # create temp Ans for betas and values to make sure order is the
    # same as Bn
    tempbetas <- c()
    tempv <- c()

    # dpsi beta1,2 for each regime dtheta
    for (j in 1:length(dimbetas)){
      # do this for each regime
      regimejdims <- dimbetas[[j]]
      dimbeta1 <- regimejdims[[1]]
      dimbeta2 <- regimejdims[[2]]
      # create the leading and trailing zero fillers for the
      # value estimating equation
      # start by updating the all
      # use the sum function to avoid errors from adding nulls
      leadbetas <- sum(leadbetas, dimbeta1, dimbeta2)
      followbetas <- sum(followbetas, -dimbeta1, -dimbeta2)
      leadv <- sum(leadv, 1)
      followv <- sum(followv, -1)
      # reset to null if j=1 or length of dimbetas
      if (j==1){
        leadbetas <- NULL
        leadv <- NULL
      } else if(j==length(dimbetas)){
        followbetas <- NULL
        followv <- NULL
      }
      # these will need to use the regime specific dataframes
      dpsibeta <- getdPsiBeta(dfValueList[[j]], q2=qFitsList[[j]]$q_fits[[2]],
                              q1=qFitsList[[j]]$q_fits[[1]], dimbeta1=dimbeta1,
                              dimbeta2=dimbeta2, leadingzeros=leadingzeros,
                              dimBn=dimBn, i=i)
      tempbetas <- rbind(tempbetas, dpsibeta)

      # now get value estimates
      dv <- getPartialValueSingleRegime(df_order=dfValueList[[j]],
                                        q2=qFitsList[[j]]$q_fits[[2]],
                                        q1=qFitsList[[j]]$q_fits[[1]],
                                        nuFits=nuFits,
                                        leadbetas=leadbetas, followbetas=followbetas,
                                        leadv=leadv, followv=followv, i=i)
      tempv <- rbind(tempv, dv)
      # update leading zeros
      leadingzeros <- leadingzeros + dimbeta1 + dimbeta2

    }
    AnTemp <- rbind(AnTemp, tempbetas, tempv)

    # there may be NaNs since pi2 is 0 for those who have not made it to that step.
    # In such cases, these should be treated as 0s.
    AnTemp[is.na(AnTemp)] <- 0
    An <- An + AnTemp
  }
  return(An/nuFits$ns)
}

#### value estimator ####

# steps are decision points
vAIPW <- function(dat, ti, steps, q_list,
                  regime_list,
                  feasibleSetsIndicator=NULL){
  df <- dat # create a copy of the data frame
  df['kappa'] <- getKappa(df=dat, ti=ti, maxEvalNumber=steps+1)
  n <- nrow(df)
  # if no one has finished do not run
  stopifnot(max(df['kappa']) == steps+1)

  # fit nu
  nuFits <- getNu(df, steps)
  # fit pi
  piFits <- getPi(df)

  qFitsList <- list()
  vList <- list()
  dfValueList <- list()

  for (regimeNumber in 1:length(regime_list)){
    regime_ind <- regime_list[[regimeNumber]]$regime_ind
    regime <- regime_list[[regimeNumber]]$regime
    # now we can turn attention to the regimes
    # add particular regime ind to dataset
    df_regimes <- cbind(df, regime_ind, regime)
    # fit regressions
    qFits <- getQfits(df=df_regimes, q_list=q_list, regime=regime,
                      feasibleSetsIndicator=feasibleSetsIndicator)
    # create the concatenated matrix to find values needed
    # do this for the regime considered

    df_value <- cbind(df, qFits[['mod_regime_vhats']],
                      qFits[['unmod_regime_vhats']],
                      piFits, regime_ind, regime)

    vTerms <- getVterms(df_value, steps, nuFits$nu)

    vi <- (vTerms['lead'] - vTerms['aug1'] - vTerms['aug2']) / nuFits$ns
    # calculate v
    vhat <- sum(vi)

    qFitsList[[regimeNumber]] <- qFits
    vList[[regimeNumber]]  <- vhat
    dfValueList[[regimeNumber]] <- vTerms
  }
  # now find the variance
  Bn <- getBn(dfValueList, qFitsList, nuFits, vList)
  # note that Bn as calculated leaves out the n-n(s) terms yet to be viewed
  # while the estimating equations for these terms are indeed 0,
  # they should be accounted for here in determining Bn
  # Bn <- (nuFits$ns / n)*Bn
  # 07052021 - the n-n(s) terms should not be used in the variance calculations.
  An <- getAn(dimBn=dim(Bn)[[1]], vList=vList,
              K=steps, dfValueList=dfValueList,
              qFitsList=qFitsList, nuFits=nuFits)

  Vn <- solve(An) %*% Bn %*% t(solve(An))

  # Values returns the estimated values for the embedded regimes
  # Vn returns the estimate of the covariance of the M-estimator vector
  return(list('Values' = vList,
              'Covariance' = Vn,
              'An' = An,
              'Bn' = Bn,
              'ns' = nuFits$ns))
}



## qstep function
#' function takes in parameters relating to a single Q step
#' the model is fit for the response passed to the function
#' predictions of the unmodified data are obtained
#' then we update the treatment to match the regime that individuals
#' should receive that is consistent
#' and obtain predictions for this data
#' We note that because of the regime consistency indicator, inconsistent
#' regimes that have modified regime here will not have their predictions used
#' as a result of the form of the estimator. 
#' We allow them to be predicted here as a means of reducing complexity 
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

## function to fit all Q functions across all K decisions
#' the feasibleSetsIndicator should be FALSE if at each decision point
#' individuals are always re-randomized, independent of response status
#' If individuals are not re-randomized (i.e. for some stage a response 
#' status prevents individuals from receiving a random treatment)
#' Then feasibleSetsIndicator should NOT be FALSE, set it to TRUE. 
getQfits <- function(df, q_list, regime, feasibleSetsIndicator=FALSE){
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
    # feasible sets indicator is a boolean that indicates if some stages have
    # only 1 treatmenat (TRUE) or multiple (FALSE)
    
    # pass previous expected outcome back
    # these will be updated by regimes in following code
    mod_regime_vhats[, k] <- mod_regime_vhats[,k+1]
    unmod_regime_vhats[, k] <- unmod_regime_vhats[,k+1]
    
    
    if (!feasibleSetsIndicator){
      # multiple options for each stage every time
      ones <- (df['kappa']>=(k+1)) # those who have finished
      new <- (df['kappa']>=k) # those who can attain a prediction from model  
      
      # now fit the model and get predictions
      vk <- qstep(qmodel = q_list[[k]],
                  data = df[ones,, drop = FALSE],
                  response = mod_regime_vhats[ones, k+1, drop = FALSE],
                  newdata = df[new,, drop = FALSE],
                  regime = regime[new,k],
                  txName = paste0('a',k) )

      # update for those who were able to receive
      mod_regime_vhats[new, k] <- vk$hats_mod
      unmod_regime_vhats[new, k] <- vk$hats_unmod
      
      # add the qfit to the outlist
      q_fits[[k]] <- vk$qfit
      
    } else{# begin when feasible sets occurs
      if (!(paste0("r", k) %in% colnames(df))){ # are we at a stage where
        # there are Rk columns? If no, (this part), then treat as normal
        # multiple options for each stage every time
        ones <- (df['kappa']>=(k+1)) # those who have finished
        new <- (df['kappa']>=k) # those who can attain a prediction from model  
        
        # now fit the model and get predictions
        vk <- qstep(qmodel = q_list[[k]],
                    data = df[ones,, drop = FALSE],
                    response = mod_regime_vhats[ones, k+1, drop = FALSE],
                    newdata = df[new,, drop = FALSE],
                    regime = regime[new,k],
                    txName = paste0('a',k) )
        
        # update for those who were able to receive
        mod_regime_vhats[new, k] <- vk$hats_mod
        unmod_regime_vhats[new, k] <- vk$hats_unmod
        
        # add the qfit to the outlist
        q_fits[[k]] <- vk$qfit
      } else{ # there are rk at this stage
        # if this is true, then we need to determine if there are multiple 
        # models to fit. This occurs for interim analyses when we have some
        # individuals who have made it to stage k, had response determined,
        # but have not finished. 
        
        
        if (length(q_list[[k]]) > 1){
          # this means that we have multiple functions to fit
          # q_fits[[k]] <- list("r0", "r1")
          
          ### start with those for rk==0 # non-responders ###
          ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
          
          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]][[paste0("r", 0)]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )
          
          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod
          
          # add the qfit to the outlist
          qr0 <- vk$qfit
          
          
          ### now do responders rk == 1 ###
          # we fit on k+1 individuals, but only use the model hats on stage k
          ones <- (df['kappa']==(k+1))*(df[paste0("r", k)] == 1) == 1
          new <- (df['kappa']==k)*(df[paste0("r", k)] == 1) == 1
          
          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]][[paste0("r", 1)]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )
          
          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod
          
          # add the qfit to the outlist
          qr1 <- vk$qfit
          
          q_fits[[k]] <- list("r0" = qr0,
                              "r1" = qr1)
          
        } else{ # this section may occur if we are estimating AIPW, no interim
          ### start with those for rk==0 # non-responders ###
          ones <- (df['kappa']>=(k+1))*(df[paste0("r", k)] != 1) == 1
          new <- (df['kappa']>=k)*(df[paste0("r", k)] != 1) == 1
          
          # now fit the model and get predictions
          # note that "response" in qstep is the outcome, not trt response
          vk <- qstep(qmodel = q_list[[k]],
                      data = df[ones,, drop = FALSE],
                      response = mod_regime_vhats[ones, k+1, drop = FALSE],
                      newdata = df[new,, drop = FALSE],
                      regime = regime[new,k],
                      txName = paste0('a',k) )
          
          # update for those who were able to receive
          mod_regime_vhats[new, k] <- vk$hats_mod
          unmod_regime_vhats[new, k] <- vk$hats_unmod
          
          # add the qfit to the outlist
          q_fits[[k]] <- vk$qfit
          
        } # end handling q_list length;
      } # end the if/else for existence of responders w/i feasible sets

      
    } # end the if/else of feasible sets, 
 
  } # end for loop
  
  colnames(mod_regime_vhats) <- paste0('q', 1:(K+1))
  colnames(unmod_regime_vhats) <- paste0(paste0('q', 1:(K+1)), '_nochange')
  return(list('q_fits' = q_fits,
              'mod_regime_vhats' = mod_regime_vhats,
              'unmod_regime_vhats' = unmod_regime_vhats))
} # end getQfits



# get coefficients from models 
getQcoefs <- function(q_all){
  coefs <- c()
  for (ell in 1:length(q_all)){
    for (k in 1:length(q_all[[ell]]$q_fits)){
      if (length(q_all[[ell]]$q_fits[[k]]) > 1){
        coefs <- c(coefs, 
                   unlist(lapply(q_all[[ell]]$q_fits[[k]], function(x) x@fitObj$coefficients)) )
      } else{
        coefs <- c(coefs,
                   q_all[[ell]]$q_fits[[k]]@fitObj$coefficients)
      } 
    }
    
  }
  return (coefs)
}

#### testing ####
# 
# ### Part I - Test the no feasible sets ###
# q2 <- modelObj::buildModelObj(model = ~ x11 + x12 + x21 + a1 + a2 + a1:a2 + x11:a1 + x12:a1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# q1 <- modelObj::buildModelObj(model = ~ x11 + x12 + a1 + x11:a1 + x12:a1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# q_list <- list(q1, q2)
# 
# df <- design1(vp=1, n=1000, seed=b*10000, s2=100, a1p=0.5, a2p=0.5, r2p = 0.4)
# 
# regimes <- list(c(0,0), c(0,1), c(1,0), c(1,1))
# regime_all <- regimelist(embRegimes = regimes,
#                          dat = df)
# t_s <- 700  # so I can do this for interim looks oor final
# 
# K <- length(q_list)
# n <- nrow(df)
# # estimate kappa
# kappa <- getKappa(df, t_s, K)
# # if no one has finished do not run
# stopifnot(max(kappa) == K+1)
# # estimate nu
# df[,'kappa'] <- kappa
# qsnf <- getQfits(df, q_list, regime_all[[1]]$regime, feasibleSetsIndicator=FALSE)
# 
# qsnf[[1]]
# 
# 
# ### Part II - Test the feasible sets, AIPW (i.e. q_list has no duplicate length) ###
# qsf <- getQfits(df, q_list, regime_all[[1]]$regime, feasibleSetsIndicator=TRUE)
# qsf[[1]] # should look different since we now are fitting only the responders
# 
# ### Part III - Test feasible sets, IAIPW (i.e. q_list has r0/r1 split) ###
# q2 <- modelObj::buildModelObj(model = ~ x11 + x12 + x21 + a1 + a2 + a1:a2 + x11:a1 + x12:a1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# q2r <- modelObj::buildModelObj(model = ~ x11 + x12 + x21 + a1 + x11:a1 + x12:a1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# q1 <- modelObj::buildModelObj(model = ~ x11 + x12 + a1 + x11:a1 + x12:a1,
#                               solver.method = 'lm',
#                               predict.method = 'predict.lm')
# q_listfr <- list(q1, list("r0" = q2,
#                         "r1" = q2r) )
# qsfr <- getQfits(df, q_listfr, regime_all[[1]]$regime, feasibleSetsIndicator=TRUE)
# qsfr[[1]]
library(dplyr)
# getPis
# only pass already subsetted data
vIPW <- function(dat, regime_list){

  getPi <- function(dat){
    # this assumes only the reduced data set is passed
    pi1_table <- dat %>% count(a1) %>% mutate(freq=n/sum(n))
    dat <- dat %>% left_join(pi1_table[,c('a1', 'freq')], by="a1")
    names(dat)[names(dat) == 'freq'] <- 'pi1'
    names(pi1_table) <- c('a1', 'n1', 'freq1')

    pi2_table <- subset(dat, kappa == 1) %>%
      count(a1, a2) %>%
      group_by(a1) %>%
      mutate(freq=n/sum(n))
    names(pi2_table) <- c('a1', 'a2', 'n2', 'pi2')
    dat <- dat %>% left_join(pi2_table[,c('a1', 'a2', 'pi2')], by=c("a1", "a2"))
    return(dat)
  }

  getV <- function(dat, ci){
    # 1/n(s) times (deltai Yi Ci) over (nu pi1 pi2)
    vhat <- sum(ci*(dat$y) / (dat$pi1*dat$pi2)) / nrow(dat)
    return(vhat)
  }

  Bn <- function(dat, vhat, regime_list, values){
    ci <- (dat$a1 == 0)*(dat$a2 == 0)
    # now find An
    # we estimated v, pis, and nu
    psi.pi1 <- (dat[,'a1'] == 1) - abs(dat[,'pi1'] - 1*(dat[,'a1'] == 0))
    # change made 6.9 based on how this is actually estimated, needs a1 as either 0 or 1
    psi.pi211 <- (dat[,'a1'] == 1)*((dat[,'a2'] == 1) - (abs(dat[,'pi2'] - (dat[,'a2'] == 0) )*(dat[,'a1'] == 1)))
    psi.pi201 <- (dat[,'a1'] == 0)*((dat[,'a2'] == 1) - (abs(dat[,'pi2'] - (dat[,'a2'] == 0) )*(dat[,'a1'] == 0)))
    # v
    psi.v <- c()
    for (regimeind in (1:length(regime_list))){
      ci <- regime_list[[regimeind]]$regime_ind[,1] * regime_list[[regimeind]]$regime_ind[,2]
      psi.v <- cbind(psi.v, (ci*(dat$y) / (dat$pi1*dat$pi2)) - values[regimeind])
    }

    psi <- as.data.frame(cbind(psi.pi1, psi.pi211, psi.pi201, psi.v))

    Bn <- (t(as.matrix(psi)) %*% as.matrix(psi)) / dim(psi)[1]
    return(Bn)
  }
  An <- function(dat, vhat, dimBn, regime_list, values){

    Antemp <- matrix(0, nrow=dimBn, ncol=dimBn)

    for (i in 1:nrow(dat)){
      Anv <- c()
      for (regimeind in (1:length(regime_list))){
        ci <- regime_list[[regimeind]]$regime_ind[i,1] * regime_list[[regimeind]]$regime_ind[i,2]
        vi <- ci*dat[i,'y'] / (dat[i, 'pi1']*dat[i,'pi2'])
        trail <- rep(0,length(values))  # this
        trail[regimeind] <- -1
        Anv <- rbind(Anv, c((-vi/dat[i,'pi1'])*(2*(dat[i,'a1'] == 1)-1),
                        (-vi/dat[i,'pi2'])*(2*(dat[i,'a1'] == 1)*(dat[i,'a2'] == 1)-1),
                        (-vi/dat[i,'pi2'])*(2*(dat[i,'a1'] == 0)*(dat[i,'a2'] == 1)-1),
                        trail) )
      }
      Antemp <- Antemp + rbind(c(-1, rep(0,dimBn-1)),
                               c(rep(0,1),-(dat[i,'a1']==1)*1, rep(0,dimBn-2)),
                               c(rep(0,2),-(dat[i,'a1']==0)*1, rep(0,dimBn-3)),
                               Anv)
    }
    return(Antemp / nrow(dat))

  }


  # find the different propensity scores
  dat <- getPi(dat)

  # estimate the different values for regimes
  values <- c()
  for (regimeind in (1:length(regime_list))){
    ci <- regime_list[[regimeind]]$regime_ind[,1] * regime_list[[regimeind]]$regime_ind[,2]
    values <- c(values, getV(dat, ci))
  }

  Bn <- Bn(dat, vhat, regime_list, values)
  An <- An(dat, vhat, dim(Bn)[[1]], regime_list, values)
  # Now find Vn
  Vn <- solve(An) %*% Bn %*% t(solve(An))
  Vn  ## the fourth diagonal is the variance of the ipw estimator vhat

  return(list('value' = values,
              'covariance' = Vn))

}


# set.seed(1)
# dat <- cr2tm(vp=1, n=1000, cont=FALSE, seed=1, s2=100, imtime = 0.5, k2time=0.5, k1time=0.5)
# dat$kappa <- (dat$t3 <= 2.5)*1  # now an indicator to be used in data
#
# dat1 <- subset(dat, kappa==1)  # find the data for the first IPW estimator
# regime_list <- regimelist(embRegimes = list(c(0,0), c(0,1), c(1,0), c(1,1)),
#                           dat = dat1)
#
# ipwres1 <- vIPW(dat1, regime_list)




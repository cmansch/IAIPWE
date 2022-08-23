#### source functions ####
library(dplyr)
library(foreach)
library(doParallel)

# source scripts 
for (file in list.files("./value_estimator", pattern=".R")){
  source(paste0("./value_estimator/", file)) # these include the data settings
}
# keep main.R out of this on the server. 
for (file in list.files("./simulation_code", pattern=".R")){
  source(paste0("./simulation_code/", file)) # these include the data settings
}
regimes <- list(list(0,c(0,0)),
                list(0,c(0,1)),
                list(0,c(1,0)),
                list(0,c(1,1)),
                list(1,c(0,0)),
                list(1,c(0,1)),
                list(1,c(1,0)),
                list(1,c(1,1)))
steps <- 2  # number of decisions K
# hypothesis testing settings
delta <- rep(22.5,8)  # the true value of regimes under the null hypothesis
H0 <- delta  
# the added part is for the non mean zero covariates 
HaltDelta <- c(10, 10, 10, 10, 0, 0, 0, 0)  # the numerator (V-H0) to be used in power calculations

t_s_list = list(500, TRUE)  # when to perform analyses, TRUE means do a final analyses
feasibleSetsIndicator=FALSE

regimelist_case_study <- function(embRegimes, dat){
  regime_list <- list()
  nRegimes <- length(embRegimes)
  n <- nrow(dat)
  K <- length(embRegimes[[1]])
  
  for(r in c(1:nRegimes)){
    regime <- matrix(NA, nrow=n, ncol=K)
    regime_ind <- matrix(NA, nrow=n, ncol=K)
    
    regimen <- embRegimes[[r]]
    for (k in 1:K){
      if (paste0("r", k) %in% colnames(dat)) {
        # then we have trt for non respond followed by trt for resp
        regime[,k] <- regimen[[k]][1]*(1-dat[,paste0("r", k)]) + 
          regimen[[k]][2]*(dat[,paste0("r", k)])
        
      } else{
        regime[,k] <- regimen[[k]]
      }
      
      regime_ind[,k] <- (dat[, paste0('a',k)] == regime[,k])*1
      
      
    }
    colnames(regime) <- paste0('regime', 1:K)
    colnames(regime_ind) <- paste0('a', 1:K, '_ind')  
    
    regime_list[[r]] <- list('regime' = regime,
                             "regime_ind" = regime_ind)
  }
  return(regime_list)
}

q2 <- modelObj::buildModelObj(model = ~ x11 + x12 + x13 + x14 + x15 +x20 + x21 +
                                a1 + a2 + a1:a2 + r2 + a1:r2,
                              solver.method = 'lm',
                              predict.method = 'predict.lm')
q1 <- modelObj::buildModelObj(model = ~ x11 + x12 + x13 + x14 + x15 +
                                a1,
                              solver.method = 'lm',
                              predict.method = 'predict.lm')
q_list <- list(q1, q2)


p1 <- modelObj::buildModelObj(model = ~ 1,
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))

p2 <- modelObj::buildModelObj(model = ~ I(a1==0):I(r2==0) -1, # + I(a1==0):I(r2==1) + I(a1==1):I(r2==0),
                              solver.method = 'glm',
                              solver.args = list(family='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))

pi_list <- list(p1, p2)

### stopping boundaries ###
#pocock bounds and OF bounds
# IAIPW
pbound <- c(2.66, 2.66)
ofbound <- c(4.20, 2.43)

# complete case AIPW
pbound_aipw <- c(2.66, 2.66)
ofbound_aipw <- c(4.30, 2.44)

# complete case IPW
pbound_ipw <- c(2.66, 2.66)
ofbound_ipw <- c(4.30, 2.44)

N <- 284  # this is based on the effective sample size for the cancer trial. 

df <- read.csv("./case_study_data.csv", row.names=1)
regime_all <- regimelist_case_study(embRegimes = regimes,
                        dat = df)

df_cc <- df[df$t3 <= 500,]
regime_all_cc <- regimelist_case_study(embRegimes = regimes,
                                    dat = df_cc)

# save the df
# write.csv(df, paste0(filedir, "/case_study_data.csv"))

#### IAIPWE ####
t1_res <- IAIPW(df=df, pi_list=pi_list, q_list=q_list, regime_all=regime_all, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=500) 
t2_res <- IAIPW(df=df, pi_list=pi_list, q_list=q_list, regime_all=regime_all, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=max(df$t3)) 

unpack_IAIPW <- function(res){
  L <- length(res$values)
  n <- res$nus$ns
  endi <- dim(res$covariance)[2]
  starti <- endi - L + 1
  cov_ <- res$covariance[starti:endi, starti:endi]
  se_ <- sqrt(diag(cov_) / n)
  
  return( list('value' = res$values,
               'variance' = cov_,
               'se' = se_,
               'n' = n))
}

t1_sum <- unpack_IAIPW(t1_res)
z1 <- (t1_sum$value - delta) / t1_sum$se  # z-statistics

t2_sum <- unpack_IAIPW(t2_res)
z2 <- (t2_sum$value - delta) / t2_sum$se  # z-statistics

z1 > pbound[1]
z1 > ofbound[1]

z2 > pbound[2]
z2 > ofbound[2]

# table data
cbind(1:8, round(t1_sum$value,1), round(t1_sum$se,1), round(z1,2), 
      round(t2_sum$value,1), round(t2_sum$se,1), round(z2,2))
# now get the z-statistics 
sum(df$t1 <= t_s_list[[1]] ) / nrow(df)
sum(df$t2 <= t_s_list[[1]] ) / nrow(df)
sum(df$t3 <= t_s_list[[1]] ) / nrow(df)

# get number of individuals
max(df$t3)
sum(df$t1 <= t_s_list[[1]] )

# get time saved if stopped early
(max(df$t3) - 500) 
(max(df$t3) - 500) / 7



#### IPWE ####
t1_resi <- IAIPW(df=df_cc, pi_list=pi_list, q_list=NULL, regime_all=regime_all_cc, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=500) 
t2_resi <- IAIPW(df=df, pi_list=pi_list, q_list=NULL, regime_all=regime_all, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=max(df$t3)) 


t1_sum <- unpack_IAIPW(t1_resi)
z1 <- (t1_sum$value - delta) / t1_sum$se  # z-statistics

t2_sum <- unpack_IAIPW(t2_resi)
z2 <- (t2_sum$value - delta) / t2_sum$se  # z-statistics

z1 > pbound[1]
z1 > ofbound[1]

z2 > pbound[2]
z2 > ofbound[2]

# table data
cbind(1:8, round(t1_sum$value,1), round(t1_sum$se,1), round(z1,2), 
      round(t2_sum$value,1), round(t2_sum$se,1), round(z2,2))
# now get the z-statistics 
sum(df$t1 <= t_s_list[[1]] ) / nrow(df)
sum(df$t2 <= t_s_list[[1]] ) / nrow(df)
sum(df$t3 <= t_s_list[[1]] ) / nrow(df)

# get number of individuals
max(df$t3)
sum(df$t1 <= t_s_list[[1]] )

# get time saved if stopped early
(max(df$t3) - 500) 
(max(df$t3) - 500) / 7


#### AIPWE ####
t1_resa <- IAIPW(df=df_cc, pi_list=pi_list, q_list=q_list, regime_all=regime_all_cc, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=500) 
t2_resa <- IAIPW(df=df, pi_list=pi_list, q_list=q_list, regime_all=regime_all, 
                feasibleSetsIndicator=feasibleSetsIndicator, t_s=max(df$t3)) 



t1_sum <- unpack_IAIPW(t1_resa)
z1 <- (t1_sum$value - delta) / t1_sum$se  # z-statistics

t2_sum <- unpack_IAIPW(t2_resa)
z2 <- (t2_sum$value - delta) / t2_sum$se  # z-statistics

z1 > pbound[1]
z1 > ofbound[1]

z2 > pbound[2]
z2 > ofbound[2]

# table data
cbind(1:8, round(t1_sum$value,1), round(t1_sum$se,1), round(z1,2), 
      round(t2_sum$value,1), round(t2_sum$se,1), round(z2,2))
# now get the z-statistics 
sum(df$t1 <= t_s_list[[1]] ) / nrow(df)
sum(df$t2 <= t_s_list[[1]] ) / nrow(df)
sum(df$t3 <= t_s_list[[1]] ) / nrow(df)

# get number of individuals
max(df$t3)
sum(df$t1 <= t_s_list[[1]] )

# get time saved if stopped early
(max(df$t3) - 500) 
(max(df$t3) - 500) / 7

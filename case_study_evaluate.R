# this is to generate a sample for the case study and estimate the outcomes

#### source functions ####

if(substr(getwd(), 2, 5) == 'home'){
  rm(list=ls())
}
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
# the local here keeps loop fast
source("./smartim/R/vipw.R", local=TRUE)
source("./smartim/R/findPocockOB.R", local=TRUE)

filedir = "/home/cmansch/elmd/sim_output/case_study"
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


#### data generation and regimes ####
design_case_study <- function(vp, n, seed=1, s2=900, a1p=0.5, a2p=0.5, r2p = 0.5){
  set.seed(seed)
  # get the parameter values associated with the value pattern given by user
  if (vp==1){
    # to vary x higher, use 0.02, 3, 3
    vpbeta <- c(-9, c(0, 0.2, 0, 10, -10), c(1), c(0), c(0, 0, 0), c(0, 0, 0))    # all equal to 22.5
  } else if (vp==2) {
    vpbeta <- c(1, c(0, 0.2, 0, 10, -10), c(1), c(0), c(-10, 0, 0), c(0, 0, 0))  # power a1=0 higher by 10
  } else if (vp==3) {
    vpbeta <- c(1, c(0, 0.2, 0, 10, -10), c(1), c(0), c(-10, -5, -2), c(10, -2, 0))
  }
  # test vps
  # a1 <- c(0,0,0,0,1,1,1,1)
  # r2 <- c(0,0,1,1,0,0,1,1)
  # a2 <- c(0,1,0,1,0,1,0,1)
  # x11 <- 152
  # x12 <- 55
  # x13 <- 0.65
  # x14 <- 0.4
  # x15 <- 0.6
  # x20 <- 22.5
  # x21 <- 0.75
  # vs <-  cbind(1, x11, x12, x13, x14, x15, x20, x21, a1, a2, a1*a2, r2, a1*r2, a1*a2*x13)  %*% vpbeta
  # mean(vs[c(1,3)])
  # mean(vs[c(1,4)])
  # mean(vs[c(2,3)])
  # mean(vs[c(2,4)])
  # mean(vs[c(5,7)])
  # mean(vs[c(5,8)])
  # mean(vs[c(6,7)])
  # mean(vs[c(6,8)])
  # 
  # initial pain should be between 5 and 10
  y1 <- as.integer(runif(n, min=50, max=101))/10 
  
  # randomly assign treatments to each individual
  a1 <- rbinom(n=n, size=1, prob=a1p)
  
  # now the design specific for covariates (mean 0)
  # note that all individuals are female (assumed) since they must have had breast cancer
  # https://academic.oup.com/biomedgerontology/article/55/11/M684/563331
  x11 <- rnorm(n=n, mean=152, sd=5) # height mean 152cm (~5ft), sd=5
  x12 <- rnorm(n=n, mean=55, sd=10)  # weight mean 55kg (~120lbs), sd=10
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6225860/
  x13 <- rbinom(n=n, size=1, p=0.65) # presence of comorbidities, maybe hypertension
  x14 <- rbinom(n=n, size=1, p=0.40) # use of pain medications (set at 40%)
  # https://www.cancer.org/content/dam/cancer-org/research/cancer-facts-and-statistics/breast-cancer-facts-and-figures/breast-cancer-facts-and-figures-2019-2020.pdf
  x15 <- rbinom(n=n, size=1, p=0.60) # if they received chemotherapy
  
  r2 <- rbinom(n=n, size=1, p=r2p)
  a2 <- rbinom(n=n, size=1, prob=a2p)
  # time 2 reduction for responders vs. non-responders
  red2 <- runif(n=n, min=30, max=40)*(r2) + runif(n=n, min=0, max=20)*(1-r2)  # reduction as percent at time2
  pain2 <- round(y1 * (1-red2/100), 1)
  # to check see summary((y1 - pain2) / y1)
  
  # adherence
  x21 <- runif(n=n, min=0.5, max=1)  # therapist adherence
  
  # errors in the reduction observed
  ei <- rnorm(n=n, mean=0, sd=sqrt(s2))
  
  # final reduction
  red3 <- cbind(1, x11, x12, x13, x14, x15, red2, x21, a1, a2, a1*a2, r2, a1*r2, a1*a2*x13) %*% 
    vpbeta  + ei
  
  # final y
  y3 <- round(y1 * (1-red3/100), 1)
  
  
  # create the times of arrival for patients
  # these are done using the fixed times given to control for variance due
  # to enrollment process.
  # We multiple the process by 1000 to have the unit be 'days'.
  t1 <- round(runif(n, 0, 1), 3)*1000
  t2 <- t1 + 56 # round(runif(n, min=35, max=57)) # move on between 5 and 8 weeks
  t3 <- t2 + 126 # round(runif(n, min=35, max=57)) + round(runif(n, min=150, max=210))
  dat <- as.data.frame(cbind(round(red3, 1), round(x11), round(x12,1), x13, x14, x15, round(red2,1), round(x21,1), a1, a2, r2, t1, t2, t3))
  names(dat)[1] <- "y"
  names(dat)[2] <- "x11"
  names(dat)[3] <- "x12"
  names(dat)[7] <- "x20"
  names(dat)[8] <- "x21"
  
  # # how many people will be through each part
  # timeana <- quantile(t3, 0.3)
  # mean(t1 > timeana)
  # mean((t1 < timeana)*(t2 > timeana))
  # mean((t1 < timeana)*(t2 < timeana)*(t3 > timeana))
  # mean((t1 < timeana)*(t2 < timeana)*(t3 < timeana))
  
  
  
  return(dat)
}

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

# seed 7 rejects pocock, not OF
# seed 17 rejects all

# df <- design_case_study(vp=3, n=N, seed=10)  # vp 3 allows the regimes to vary from 19.5 to 37.5
# # save the df
# write.csv(df, paste0(filedir, "/case_study_data.csv"))

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





#### look at ses

t1_sum <- unpack_IAIPW(t1_res)
t1_suma <- unpack_IAIPW(t1_resa)
t1_sumi <- unpack_IAIPW(t1_resi)
cbind(t1_sum$se, t1_suma$se, t1_sumi$se)  # not hugely more efficient?


t2_sum <- unpack_IAIPW(t2_res)
t2_suma <- unpack_IAIPW(t2_resa)
t2_sumi <- unpack_IAIPW(t2_resi)
cbind(t2_sum$se, t2_suma$se, t2_sumi$se)  # not hugely more efficient?

z1 <- (t1_sum$value - delta) / t1_sum$se  # z-statistics
z2 <- (t2_sum$value - delta) / t2_sum$se  # z-statistics
z1a <- (t1_suma$value - delta) / t1_suma$se  # z-statistics
z2a <- (t2_suma$value - delta) / t2_suma$se  # z-statistics
z1i <- (t2_sumi$value - delta) / t1_sumi$se  # z-statistics
z2i <- (t2_sumi$value - delta) / t2_sumi$se  # z-statistics

out_table <- cbind(1:8, 
      round(t1_sumi$value,1) / 10, round(t1_sumi$se,2), round(z1i, 2),
      round(t1_suma$value,1) / 10, round(t1_suma$se,2), round(z1a, 2),
      round(t1_sum$value,1) / 10, round(t1_sum$se,2), round(z1, 2),
      round(t2_sumi$value,1) / 10, round(t2_sumi$se,2), round(z2i, 2),
      round(t2_suma$value,1) / 10, round(t2_suma$se,2), round(z2a, 2),
      round(t2_sum$value,1) / 10, round(t2_sum$se,2), round(z2, 2)
)
      

out_table <- cbind(1:8, 
      round(t1_sumi$value,1) / 10, round(t1_sumi$se,2),
      round(t1_suma$value,1) / 10, round(t1_suma$se,2),
      round(t1_sum$value,1) / 10, round(t1_sum$se,2),
      round(t2_sumi$value,1) / 10, round(t2_sumi$se,2),
      round(t2_suma$value,1) / 10, round(t2_suma$se,2),
      round(t2_sum$value,1) / 10, round(t2_sum$se,2)
)
out_table
apply(out_table, 1, paste, collapse=" & ")


out_table <- cbind(1:8,
      pbound[1], ofbound[1],
      round(t1_sum$value,1) / 10, paste0("(", round(t1_sum$se,2), ")"), round(z1, 2),
      pbound[2], ofbound[2],
      round(t2_sum$value,1) / 10, paste0("(", round(t2_sum$se,2), ")"), round(z2, 2)
)
apply(out_table, 1, paste, collapse=" & ")


# get number of individuals
max(df$t3)
sum(df$t1 <= t_s_list[[1]] ) / nrow(df)
sum(df$t2 <= t_s_list[[1]] ) / nrow(df)
sum(df$t3 <= t_s_list[[1]] ) / nrow(df)

# get time saved if stopped early
(max(df$t3) - 500) 
(max(df$t3) - 500) / 7


out_tablea <- cbind(1:8,
                   pbound_aipw[1], ofbound_aipw[1],
                   round(t1_suma$value,1) / 10, paste0("(", round(t1_suma$se,2), ")"), round(z1a, 2),
                   pbound_aipw[2], ofbound_aipw[2],
                   round(t2_suma$value,1) / 10, paste0("(", round(t2_suma$se,2), ")"), round(z2a, 2)
)
apply(out_tablea, 1, paste, collapse=" & ")



out_tablei <- cbind(1:8,
                    pbound_ipw[1], ofbound_ipw[1],
                    round(t1_sumi$value,1) / 10, paste0("(", round(t1_sumi$se,2), ")"), round(z1i, 2),
                    pbound_ipw[2], ofbound_ipw[2],
                    round(t2_sumi$value,1) / 10, paste0("(", round(t2_sumi$se,2), ")"), round(z2i, 2)
)
apply(out_tablei, 1, paste, collapse=" & ")

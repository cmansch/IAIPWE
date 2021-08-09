#' this function finds the required sample sizes and n1, n2 splits for the
#' case when interim analyses occur for only completed cohorts
#' It uses the scripts vp1, vp2, and vp3 in addition to the findPocockOB
#' script to find the boundaries and the sample sizes required for the
#' simulations.
#' Note, it also finds the case in which there is no interim analysis so that
#' there is a base to compare against.

library(smartim)

# vp1cor <- as.matrix(read.csv(paste0(getwd(), "/sim_output/vp1_Zcor.csv"), row.names=1))
# vp1var <- as.matrix(read.csv(paste0(getwd(), "/sim_output/vp1_var.csv"), row.names=1))

enrolls <- list()

enrolls[[1]] <- c(0.5,0.6,0.7)
enrolls[[2]] <- c(0.5,0.7,0.8)
enrolls[[3]] <- c(0.5,0.8,0.9)
enrolls[[4]] <- c(0.4,0.5,0.6)

for (enroll in enrolls){
  im <- enroll[1]
  k2 <- enroll[2]
  k1 <- enroll[3]
  
  vp <- 1
  setting <- paste('vp', vp, 'im', im, 'k2', k2, 'k1', k1)
  vp1cor <- as.matrix(read.csv(paste0(getwd(), "/sim_output/", setting, "_Zcor.csv"), row.names=1))
  vp1var <- as.matrix(read.csv(paste0(getwd(), "/sim_output/", setting, "_var.csv"), row.names=1))
  
  n2split <- max(im, k2, k1)
  
  # find the boundaries
  vp1P <- findPocockBoundary(corr=vp1cor, alpha=0.05, startvalue = 2)
  vp1OF <- findOFBoundary(corr=vp1cor, alpha=0.05, startvalue = 2)
  
  vp1Fixed <- findPocockBoundary(corr=vp1cor[(dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1]),
                                             (dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1])],
                                 alpha=0.05, startvalue = 2)
  # find the fixed sample size
  vp2FixedN <- findFixedSample(vp1cor[((dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1])), 
                                      ((dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1]))], 
                               pow=0.8,
                               cov=vp1var[(dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1])],
                               vp=c(0,0,0,2), bound=vp1Fixed$bound,
                               lambda=10, N0=50)
  
  vp3FixedN <- findFixedSample(vp1cor[((dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1])), 
                                      ((dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1]))], 
                               pow=0.8,
                               cov=vp1var[(dim(vp1cor)[1]/2+1) :(dim(vp1cor)[1])],
                               vp=c(0,2,2,0), bound=vp1Fixed$bound,
                               lambda=10, N0=50)
  
  # find the Pocock and OF sample sizes required to attain 0.8 power
  # under vp2
  
  vp2PocockN <- findPocockSample(vp1cor, 0.8, vp1var, vp=c(0,0,0,2,0,0,0,2),
                                 bound=vp1P$bound, interim=n2split)
  
  vp2OFN <- findOFSample(vp1cor, 0.8, vp1var, vp=c(0,0,0,2,0,0,0,2),
                         bound=vp1OF$bound, interim=n2split)
  
  # under vp3
  
  vp3PocockN <- findPocockSample(vp1cor, 0.8, vp1var, vp=c(0,2,2,0,0,2,2,0),
                                 bound=vp1P$bound, interim=n2split)
  
  vp3OFN <- findOFSample(vp1cor, 0.8, vp1var, vp=c(0,2,2,0,0,2,2,0),
                         bound=vp1OF$bound, interim=n2split)
  
  output <- list('Fixed' = list('zcutoff' = vp1Fixed,
                                'vp2N' = vp2FixedN,
                                'vp3N' = vp3FixedN),
                 'Pocock' = list('zcutoff' = vp1P,
                                 'vp2N' = vp2PocockN,
                                 'vp3N' = vp3PocockN),
                 'OF' = list('zcutoff' = vp1OF,
                             'vp2N' = vp2OFN,
                             'vp3N' = vp3OFN) )
  out_location <-  paste0(getwd(), "/sim_output/", setting, "boundary_and_sample_sizes.rdata")
  rlist::list.save(output, file = out_location)
  
  print(output)  
}

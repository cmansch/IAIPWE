# Simulate data under a specified enrollment process.
# design 1
design1 <- function(vp, n, seed=1, s2=100){
  set.seed(seed)
  # get the parameter values associated with the value pattern given by user
  if (vp==1){
    vpbeta <- c(10, 0, 0, 0)
  } else if (vp==2) {
    vpbeta <- c(10, 0, 0, 2)
  } else if (vp==3) {
    vpbeta <- c(10, 2, 2, -4)
  }
  # randomly assign treatments to each individual
  a1 <- rbinom(n=n, size=1, prob=0.5)
  a2 <- rbinom(n=n, size=1, prob=0.5)
  # create the error for individuals
  ei <- rnorm(n=n, mean=0, sd=sqrt(s2))
  # create the observed outcome based
  designMatrix <- cbind(1, a1, a2, a1*a2)
  y <- (designMatrix %*% vpbeta) + ei

  # create the times of arrival for patients
  # these are done using the fixed times given to control for variance due
  # to enrollment process.
  # We multiple the process by 1000 to have the unit be 'days'.
  t1 <- round(runif(n, 0, 1), 3)*1000
  t2 <- t1 + 0.1*1000 + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  t3 <- t2 + 0.1*1000 + round(runif(n, -1, 1), 1)*10 # 0.1 means that 10% move on each 100 day +- 10
  dat <- as.data.frame(cbind(a1, a2, y, t1, t2, t3))
  names(dat)[3] <- "y"

  return(dat)
}

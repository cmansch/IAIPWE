#' This portion of the package is designed to make samples based on the
#' corresponding paper, data considered therein, and their value patterns.
#' See doi:
#' These functions create the designs referenced in the paper.

#' CR2
#' @param vp value pattern number 1, 2, or 3 from the paper
#' @param n the total sample size
#' @param cont boolean indicating whether to include a control arm (Y=T)
#' @param seed int, seed for randomization
#' @param s2 float, the variance of the patient outcome
#' @export
cr2 <- function(vp, n, cont=FALSE, seed=1, s2=100){
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
  # we include the addition of the maximum t1, t2
  # this design has rolling admission with some amount of randomness in entry;
  # individuals arrive at different time intervals
  # follow ups are 3-4 weeks. Total time of trial is <2 months.
  t1 <- (1:n) + sample(1:10, n, replace=TRUE)
  # t2 <- t1 + 60 + sample(-5:30, n, replace=TRUE)
  # t3 <- t2 + 60 + sample(-5:30, n, replace=TRUE)
  t2 <- t1 + sample(20:30, n, replace=TRUE)
  t3 <- t2 + sample(20:30, n, replace=TRUE)
  dat <- as.data.frame(cbind(a1, a2, y, t1, t2, t3))
  names(dat)[3] <- "y"

  return(dat)
}

#imtime is the proportion of finishers when the data is analyzed
#k1time is the number of people enrolled in the trial when the final is reached
#k2time is the number of pepole who have received the first and second treatment when final time reached
# as such imtime < k2time < k1time
cr2tm <- function(vp, n, cont=FALSE, seed=1, s2=100, imtime = 0.5, k2time=0.6, k1time=0.75){
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
  t1 <- 0 + ((1:n)/n > k1time)*1  + ((1:n)/n > k2time)*1 + ((1:n)/n > imtime)*1
  t2 <- t1 + 1
  t3 <- t2 + 1
  dat <- as.data.frame(cbind(a1, a2, y, t1, t2, t3))
  names(dat)[3] <- "y"

  return(dat)
}

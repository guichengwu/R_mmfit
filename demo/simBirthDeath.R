# First half of this file is the mathematical method. The second half, is the Simulation code.
# Mathematical Method developed in Problem 1
calcCov <- function(piVector, alphaVector) {
  sum1 <- 0
  for (i in 1:length(piVector)) {
    sum1 <- sum1 + piVector[i]*(i-1)^2
  }
  
  sum2 <- 0
  for (i in 1:length(piVector)) {
    sum2 <- sum2 + piVector[i]*(i-1)
  }
  sum2 <- sum2^2
  
  sum3 <- 0 
  for (i in 1:length(piVector)) {
    sum3 <- sum3 + piVector[i]*(i-1)*(alphaVector[i]-(1-alphaVector[i]))
  }
  
  sum4 <- 0
  for (i in 1:length(piVector)) {
    sum4 <- sum4 + piVector[i]*(i-1)
  }
  
  sum5 <- 0
  for (i in 1:length(piVector)) {
    sum5 <- sum5 + piVector[i]*(alphaVector[i]-(1-alphaVector[i]))
  }
  
  return(sum1 - sum2 + sum3 - (sum4 * sum5))
}

n <- 10000  # number of possible of states
m <- 500000 # number of time intervals
t <- 5 # number of experiments

alpha <- matrix(rep(0,n*t), ncol=n)
current=matrix(rep(0,m*t),ncol=m)
current_shifted=matrix(rep(0,m*t),ncol=m)
cumulative_states <- matrix(rep(0,n*t), ncol=n)
pi=matrix(rep(0,n*t),ncol=n)

for (i in 1:t) {
  for (j in 1:n) {
    # Considering that alphas (and betas) are uniformly distributed
    alpha[i,j] <- runif(1,0,1) 
  }
}

for (i in 1:t) {
  current_state <- 0
  for (j in 1:m) {
    current[i,j]=current_state
    cumulative_states[i, current_state + 1] <- cumulative_states[i, current_state + 1] + 1
    jump <- runif(1)
    if (jump < alpha[i, current_state + 1]) {
      if (current_state != n-1) {
        current_state <- current_state + 1
      } 
    } else {
      if (current_state != 0) {
        current_state <- current_state - 1
      }
    }
    if (j > 1){
      current_shifted[i,j]=current[i,j-1] # Track previous state
    }
  }
}

for (i in 1:t){
  for (j in 1:n){
    pi[i,j]=(cumulative_states[i,j])/m
  }
}

averageCumulativeStates <- colMeans(cumulative_states)
pis <- averageCumulativeStates/m
meanAlphas <- colMeans(alpha)
meanBetas <- colMeans(1-alpha)

for (i in 1:t) {
  simulatedCov <- cov(current[i,],current_shifted[i,])
  calculatedCov <- calcCov(pi[i,],alpha[i,])
  cat("Experiment", i,">> Simulation Cov = ", simulatedCov," | Mathematical Cov = ", calculatedCov, "| Difference = ", abs(1-calculatedCov/simulatedCov)*100,"%\n")
}



###########################################################################################################
# begin: test function

#test Poisson esitimated result
testmmfitPoisson = function() {
  poissonData <- poisson_sim(10000, 5)
  start <- c(2)
  #result <- poisson_func(start, poissonData)
  testResult <- mmfit(poisson_func, poissonData, start, 'poisson')
  return(testResult)
}

#test power law estimated result
testmmfitPowerlaw = function() {
  powerlawData = power_law_sim(10000, 1, 10)
  start <- 0.5
  #result <- power_law_func(start, powerlawData)
  testResult <- mmfit(power_law_func, powerlawData, start, 'powerlaw')
  return(testResult)
}

#test gamma estimated result
testmmfitGamma = function () {
  gammaData <- gamma_sim(10000, 7.5, 4)
  start <- c(7, 2)
  testResult <- mmfit(gamma_func, gammaData, start, 'gamma')
  
  return(testResult)
}

#test beta estimated result
testmmfitBeta = function() {
  betaData <- beta_sim(10000, 7.5, 4)
  start <- c(5, 2)
  testResult <- mmfit(beta_func, betaData, start, 'beta')
  
  return(testResult)
}

# test poisson estimated result
testmmfitPoissonMixture = function() {
  poissonMixtureData <- poisson_mixture_sim(10000, 6, 20, 0.3)
  start <- c(4, 17, 0.2)
  testResult <- mmfit(poisson_mixture_func, poissonMixtureData, start, 'poissonMixture')
  
  return(testResult)
}

# test plot estimated result
testmmfitExpMixture = function() {
  expMixtureData <- exp_mixture_sim(10000, 3, 10, 0.3)
  start <- c(2, 8, 0.2)
  testResult <- mmfit(exp_mixture_func, expMixtureData, start, 'expMixture')
  
  return(testResult)
}


#test plot poisson 
testPlotPoisson = function() {
  poissonData <- poisson_sim(10000, 8)
  start <- c(1)
  plot.mmfit(poisson_func, poissonData, start, 'poisson')
}

# test power law 
testPlotPowerlaw = function() {
  powerlawData <- power_law_sim(1000, 1, 3)
  start <- c(2)
  plot.mmfit(power_law_func, powerlawData, start, 'powerlaw')
}

#test plot gamma 
testPlotGamma = function() {
  gammaData <- gamma_sim(1000, 2, 0.5)
  start <- c(1, 0.2)
  
  plot.mmfit(gamma_func, gammaData, start, 'gamma')
}

testPlotUserFunc1 = function() {
  poissonData <- poisson_sim(10000, 8)
  start <- c(lambda = 1)
  plot.mmfit(poisson_func, poissonData, start, '', poisson_density)
}

testPlotUserFun2 = function() {
  gammaData <- gamma_sim(1000, 2, 0.5)
  start <- c('shape' = 1, 'rate' = 0.2)
  
  plot.mmfit(gamma_func, gammaData, start, 'user', gamma_density)
}

#test plot beta 
testPlotBeta = function() {
  betaData <- beta_sim(1000, 4, 3)
  start <- c(0.7, 1)
  
  plot.mmfit(beta_func, betaData, start, 'beta')
  #mmfit(beta_func, betaData, start, 'beta')
}

# test plot poisson mixture
testPlotPoisMixture = function() {
  mixtureData <- poisson_mixture_sim(100000, 6, 20, 0.3)
  start <- c(5, 19, 0.2)
  
  plot.mmfit(poisson_mixture_func, mixtureData, start, 'poissonMixture')
}

# test plot exponential mxiture
testPlotExpMixture = function() {
  expMixtureData <- exp_mixture_sim(100000, 5, 8, 0.3)
  start <- c(6, 7, 0.5)
  
  plot.mmfit(exp_mixture_func, expMixtureData, start, 'expMixture')
}

testPrintPoisson = function() {
  poissonData <- poisson_sim(10000, 8)
  start <- c(1)
  print.mmfit(poisson_func, poissonData, start, 'poisson')
}

testPrintPowerlaw = function() {
  powerlawData <- power_law_sim(10000, 1, 4)
  start <- c(k = 3)
  print.mmfit(power_law_func, powerlawData, start, 'powerlaw')
}

testPrintGamma = function() {
  gammaData <- gamma_sim(10000, 7.5, 4)
  start <- c(alpha = 5, beta = 2)
  
  print.mmfit(gamma_func, gammaData, start, 'gamma')
}

testPrintBeta = function() {
  betaData <- beta_sim(10000, 2, 4)
  start <- c(alpha = 0.5, beta = 1)
  
  print.mmfit(beta_func, betaData, start, 'beta')
}

testPrintPoisMixture = function() {
  mixtureData <- poisson_mixture_sim(100000, 6, 20, 0.3)
  start <- c(lambda1 = 5, lambda2 = 19, r = 0.2)
  
  print.mmfit(poisson_mixture_func, mixtureData, start, 'poissonMixture')
}

testPrintExpMixture = function() {
  expMixtureData <- exp_mixture_sim(100000, 0.1, 0.1, 0.3)
  start <- c(lambda1 = 0.1, lambda2 = 0.3, r = 0.4)
  
  print.mmfit(exp_mixture_func, expMixtureData, start, 'expMixture')
}

testSummaryPoisson = function() {
  poissonData <- poisson_sim(10000, 8)
  start <- c(1)
  summary.mmfit(poisson_func, poissonData, start, 'poisson')
}
# end: test function
###########################################################################################################
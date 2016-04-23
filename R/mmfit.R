###########################################################################################################
# begin: include library

library(poweRlaw)
library(ggplot2)
library(gmm)
library(plotly)

# end: include library
###########################################################################################################

###########################################################################################################
#begin: main mmfit function
#x: a vector or matrix of observations, each row contains the data for one observation
#g: a function specify the moment equations, the argument g in gmm
#start is a vector of initial guess values
mmfit = function(g, x, start, distributionType, user_density_func = "") {
  result <- 0;
  if (distributionType == 'poisson') {
    result <- summary(gmm(g, x, start, lower = 0, upper = 10, method = "Brent"))
  } 
  else if (distributionType == 'powerlaw') {
    result <- summary(gmm(g, x, start, lower = 0, upper = 10, method = "Brent"))
  }
  else if (distributionType == 'gamma' || distributionType == 'beta' || 
           distributionType == 'expMixture' || distributionType == 'poissonMixture') {
    result <- summary(gmm(g, x, start))
  }
  else {
    if (length(start) == 1) {
      result <- summary(gmm(g, x, start, lower = 0, upper = 10, method = "Brent"))
    } else {
      result <- summary(gmm(g, x, start))
    } 
  }
  
  thetahat <- result$coefficients[,1]     # thetahat, the vector of estimated parameters
  thetahatses <- result$coefficients[, 2]  #thetahatses, the vector of standard errors
  #denscomp <- list(x, thetahat) # denscomp, a graph object that compares the parametric and nonparametric density estimates
  denscomp <- 0
  cdfband <- 0   #cdfband, a graph object that draws the empirical cdf and an enclosing Kolmogorov-Smirnov confidence band
  
  #result <- mmfit(g, x, start, distributionType)
  #denscomp <- result$denscomp
  sampleDataSize <- length(x)
  sampleData <- data.frame(x)
  
  if (distributionType == 'poisson') {
    thetahat <- c(lambda = thetahat)
    thetahatses <- c(SE.lambda = thetahatses)
    denscomp <- ggplot(sampleData, aes(x)) +
      geom_histogram(aes(y=..density..),
                     binwidth = 1, 
                     colour = "black", fill = "white") +
      stat_function(fun = poisson_density, 
                    args = list(lambda = thetahat[1]),
                    colour = "blue", geom = "line") +
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
  }
  
  else if (distributionType == 'powerlaw') {
    thetahat <- c(gamma = thetahat)
    thetahatses <- c(SE.gamma = thetahatses)
    k=1
    currentsum=0
    while(TRUE){
      diff=k^(-thetahat)
      currentsum=currentsum + diff
      k=k+1
      if (diff < 0.01){
        break
      }
      c = 1/ currentsum
    }
    denscomp <- ggplot(sampleData, aes(x)) + 
      geom_histogram(aes(y=..density..),
                     binwidth=1,
                     colour = "black", fill="white") +
      stat_function(fun = powerlaw_density,
                    args = list(theta = thetahat[1], c),
                    colour = "blue", geom = "line") + 
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
  }
  
  else if (distributionType == 'gamma') {
    thetahat <- c(shape = thetahat[[1]], rate = thetahat[[2]])
    thetahatses <- c(SE.shape = thetahatses[[1]], SE.rate = thetahatses[[2]])
    denscomp <- ggplot(sampleData, aes(x)) +
      geom_histogram(aes(y=..density..),      
                     bins = 20,
                     colour="black", fill="white") +
      stat_function(fun = gamma_density, n= sampleDataSize, 
                    args = list(shape =thetahat[1], rate = thetahat[2]),
                    colour = "blue", geom = "line") + 
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
    
  }
  
  else if (distributionType == 'beta') {
    thetahat <- c(alpha = thetahat[[1]], beta = thetahat[[2]])
    thetahatses <- c(SE.alpha = thetahatses[[1]], SE.beta = thetahatses[[2]])
    denscomp <- ggplot(sampleData, aes(x)) +
      geom_histogram(aes(y=..density..),
                     bins=20,
                     colour = "black", fill="white") +
      stat_function(fun = dbeta,
                    args = list(shape1 = thetahat[1], shape2 = thetahat[2]),
                    colour = "blue", geom = "line") + 
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
  }
  
  else if (distributionType == 'expMixture') {
    thetahat <- c(lambda1 = thetahat[[1]], lambda2 = thetahat[[2]], r = thetahat[[3]])
    thetahatses <- c(SE.lambda1 = thetahatses[[1]], SE.lambda2 = thetahatses[[2]], SE.r = thetahatses[[3]])
    denscomp <-ggplot(sampleData, aes(x)) +
      geom_histogram(aes(y=..density..),      
                     bins = 20,
                     colour="black", fill="white") +
      stat_function(fun = expMixture_density, n= sampleDataSize, 
                    args = list(lambda1 = thetahat[1], lambda2 = thetahat[2], r = thetahat[3]),
                    colour = "blue", geom = "line") + 
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
  }
  
  else if (distributionType == 'poissonMixture') {
    thetahat <- c(lambda1 = thetahat[[1]], lambda2 = thetahat[[2]], r = thetahat[[3]])
    thetahatses <- c(SE.lambda1 = thetahatses[[1]], SE.lambda2 = thetahatses[[2]], SE.r = thetahatses[[3]])
    denscomp <- ggplot(sampleData, aes(x)) +
      geom_histogram(aes(y=..density..),      
                     bins=20,
                     colour="black", fill="white") +
      stat_function(fun = poisMixture_density, n= sampleDataSize, 
                    args = list(lambda1 = thetahat[1], lambda2 = thetahat[2], r = thetahat[3]),
                    colour = "blue", geom = "line") +       
      ggtitle("Parametric and nonparametric density estimate") +
      theme(plot.title=element_text(face="bold", size=10))
  }else{
    startLength <- length(start)
    startNames <- names(start)
    for (i in 1:startLength) {
      names(thetahat)[i] <- startNames[i]
      names(thetahatses)[i] <- paste('SE.', startNames[i], sep="")
    }
    
    
    # thetahat <- c(names(start[1]) = thetahat[[1]], lambda2 = thetahat[[2]], r = thetahat[[3]])
    # thetahatses <- c(SE.lambda1 = thetahatses[[1]], SE.lambda2 = thetahatses[[2]], SE.r = thetahatses[[3]])
    if (startLength == 1) {
      denscomp <- ggplot(sampleData, aes(x)) +
        geom_histogram(aes(y=..density..),      
                       bins=20,
                       colour="black", fill="white") +
        stat_function(fun = user_density_func, n= sampleDataSize, 
                      #args = list(shape = thetahat[1], rate = thetahat[2]),
                      args = list(thetahat[1]),
                      colour = "blue", geom = "line") +       
        ggtitle("Parametric and nonparametric density estimate") +
        theme(plot.title=element_text(face="bold", size=10))
    }
    else if (startLength == 2) {
      denscomp <- ggplot(sampleData, aes(x)) +
        geom_histogram(aes(y=..density..),      
                       bins=20,
                       colour="black", fill="white") +
        stat_function(fun = user_density_func, n= sampleDataSize, 
                      #args = list(shape = thetahat[1], rate = thetahat[2]),
                      args = list(thetahat[1], thetahat[2]),
                      colour = "blue", geom = "line") +       
        ggtitle("Parametric and nonparametric density estimate") +
        theme(plot.title=element_text(face="bold", size=10))
    } 
    else if (startLength == 3) {
      denscomp <- ggplot(sampleData, aes(x)) +
        geom_histogram(aes(y=..density..),      
                       bins=20,
                       colour="black", fill="white") +
        stat_function(fun = user_density_func, n= sampleDataSize, 
                      #args = list(shape = thetahat[1], rate = thetahat[2]),
                      args = list(thetahat[1], thetahat[2], thetahat[3]),
                      colour = "blue", geom = "line") +       
        ggtitle("Parametric and nonparametric density estimate") +
        theme(plot.title=element_text(face="bold", size=10))
    }
    
  }
  
  x = sort(x)
  ecdfValue = ecdf(x)
  
  if (distributionType == 'poisson' || distributionType == 'poissonMixture' || distributionType == 'powerlaw'){
    xValue=c(rep(0,(2*(max(x)-min(x)))));
    yValue=c(rep(0,(2*(max(x)-min(x)))));
    upperValue <- rep(0, length(yValue));
    lowerValue <- rep(0, length(yValue));
    
    for (i in 1:max(x)-min(x)){
      if (i==1){
        xValue[2*i]=min(x)+i
        xValue[2*i+1]=min(x)+i
        xValue[i]=min(x)
        
      }else if(i<max(x)-min(x)){
        xValue[2*i]=min(x)+i
        xValue[2*i+1]=min(x)+i
      }else{
        xValue[2*i]=max(x)
      }
    }
    for (j in 1: (max(x)-min(x))){
      yValue[2*j]=ecdfValue(xValue[2*j])
      yValue[2*j-1]=ecdfValue(xValue[2*j])
    }
    }else{
      xValue=x
      yValue <- rep(0, sampleDataSize)
  upperValue <- rep(0, sampleDataSize)
  lowerValue <- rep(0, sampleDataSize)
  
  yValue=ecdfValue(x);
    }
  upperValue <- yValue + 1.358 * ((sampleDataSize)^(-0.5));
  lowerValue <- yValue - 1.358 * ((sampleDataSize)^(-0.5));
  
  cdfDataFrame <- data.frame(x = xValue, y = yValue, yupper = upperValue, ylower = lowerValue)
  cdfband <- ggplot(cdfDataFrame) + geom_line(aes(x=xValue, y = yValue)) + 
  geom_ribbon(aes(x=xValue, ymin=ylower, ymax=yupper,fill = "band"), alpha = 0.3) +
    scale_color_manual("", values = "blue") +
    scale_fill_manual("" , values = "grey12") +
    ggtitle("Empirical cdf with K-S confidence band") +
    theme(plot.title=element_text(face="bold", size=10)) +
    labs(y= "ecdf") 
  
  mmf <- list(thetahat = thetahat, thetahatses = thetahatses, denscomp = denscomp, cdfband = cdfband)
  
  class(mmf) <- "mmf"
  
  return(mmf)
}
# end: mmfit function
###########################################################################################################

###########################################################################################################
#begin : simulate data
# simulate n samples of poisson data
poisson_sim = function(n, lambda) {
  return(rpois(n, lambda))
}

# simulate n samples of power law data
power_law_sim = function(n,xmin,k){
  p=rpldis(n,xmin,k)
  return(p)
}

# simulate n samples of gamma data
gamma_sim = function(n, shape, rate) {
  return(rgamma(n, shape, rate))
}

# simulate n samples of beta data
beta_sim = function(n, shape1, shape2) {
  return(rbeta(n, shape1, shape2))
}

# simulate one sample of mixture poisson data
poisson2 = function(lambda1,lambda2,r){
  m <- sample(1:10, 1)
  
  if (m <= r*10) {
    data <- rpois(1, lambda1)
  } else {
    data <- rpois(1, lambda2)
  }
  
  return(data)
}

#simulate n samples of mixture poisson data
poisson_mixture_sim = function(n,lambda1,lambda2,r){
  data <- c(rep(0, n))
  for (i in 1:n) {
    data[i] <- poisson2(lambda1, lambda2, r)
  }
  return(data)
}

#simulate one sample of mixture exponential data
exponential = function(lambda1, lambda2, r) {
  m <- sample(1:10, 1)
  if (m <= r*10) {
    data <- rexp(1, lambda1)
  } else {
    data <- rexp(1, lambda2)
  }
  return(data)
}

# simulate n samples of mixture exponential data
exp_mixture_sim = function(n, lambda1, lambda2, r) {
  data <- c(rep(0, n))
  for (i in 1:n) {
    data[i] <- exponential(lambda1, lambda2, r)
  }
  return(data)
}
# end: simulate data
###########################################################################################################

###########################################################################################################
# begin moment function
# poisson estimate function
poisson_func = function(theta,x){
  #   lambda=theta
  #   mean=lambda
  #   f=mean(x)
  #   return(f)
  cbind(theta - x)
}

# power law estimate function
power_law_func = function(theta,x){
  mean = (theta - 1) / (theta - 2)
  cbind(mean - x)
  
}

# gamma moment function
gamma_func = function(theta,x){
  t1=theta[1]
  t2=theta[2]
  meanb=t1/t2
  m1=meanb-x
  m2=t1/(t2^2)-(x-meanb)^2
  f=cbind(m1,m2)
  return(f)
}

# beta moment function
beta_func = function(theta,x){
  t1=theta[1]
  t2=theta[2]
  t12=t1+t2
  meanb=t1/t12
  m1 = meanb-x
  m2=t1*t2/(t12^2 * (t12+1)) - (x-meanb)^2
  f=cbind(m1,m2)
  return(f)
}

# poisson mixture moment function
poisson_mixture_func = function(theta,x){
  lambda1=theta[1]
  lambda2=theta[2]
  r=theta[3]
  mean = r*lambda1 + (1-r)*lambda2
  m1=mean-x
  m2 = r*(lambda1^2 + lambda1) + (1-r)*(lambda2^2 + lambda2) - x^2
  m3= r*(lambda1^3+3*lambda1^2+ lambda1) + (1-r)*(lambda2^3 + 3*lambda2^2 + lambda2) - x^3
  f=cbind(m1,m2,m3)
  return(f)
}

# expential mixture moment function
exp_mixture_func = function(theta,x){
  lambda1=theta[1]
  lambda2=theta[2]
  r=theta[3]
  mean = r/lambda1 + (1-r)/lambda2
  m1=mean-x
  m2 = 2*r / (lambda1^2) + 2*(1-r) / (lambda2^2) - x^2
  m3= 6*r / (lambda1^3) + 6*(1-r) / (lambda2^3) -x^3
  f=cbind(m1,m2,m3)
  return(f)
}
# end : moment function
###########################################################################################################

###########################################################################################################
#begin: density function
#poisson density function
poisson_density = function(x, lambda) {
  #dpois(x, 1.5)
  (lambda^x)*(exp(-lambda)) / factorial(x)
}

# power law density function
powerlaw_density = function(x, theta, c) {
  (c)*(x^(-theta))
}

# gamma density function
gamma_density = function(x, shape, rate) {
  ((rate^shape) * (x^(shape-1))*exp((-rate*x))) / factorial(shape - 1)
}

# beta density function
beta_density = function(x, shape1, shape2) {
  dbeta(x, shape1, shape2)
}

# expoential mixture density function
expMixture_density = function(x, lambda1, lambda2, r) {
  r * lambda1 * (exp((-lambda1*x))) + (1 - r) * lambda2 * (exp((-lambda2*x)))
}

#poisson mixture density function
poisMixture_density = function(x, lambda1, lambda2, r) {
  r * (lambda1^x)*(exp(-lambda1)) / factorial(x) + (1-r)*(lambda2^x)*(exp(-lambda2)) / factorial(x) 
}


#end: density function
###########################################################################################################

###########################################################################################################
#begin : print.mmfit(), summary.mmfit(), plot.mmfit()
print = function(x, ...) {
  UseMethod("print", x)
}

plot = function(x, ...) {
  UseMethod("plot", x)
}

summary = function(x, ...) {
  UseMethod("summary", x)
}

print.mmfit <- function(g, x, start, distributionType, user_density_func) {
  result <- mmfit(g, x, start, distributionType, user_density_func)
  thetahat <- c(result$thetahat)
  thetahatses <- c(result$thetahatses)
  print(thetahat)
  print(thetahatses)
}

summary.mmfit = function(g, x, start, distributionType, user_density_func) {
  print.mmfit(g, x, start, distributionType, user_density_func)
  plot.mmfit(g, x, start, distributionType, user_density_func)
}

plot.mmfit = function(g, x, start, distributionType, user_density_func) {
  require(gridExtra)
  result <- mmfit(g, x, start, distributionType, user_density_func)
  denscomp <- result$denscomp
  cdfband <- result$cdfband
  
  grid.arrange(denscomp, cdfband, ncol = 2)
  #ggplotly()
}
#end: print.mmfit(), summary.mmfit(), plot.mmfit()
###########################################################################################################



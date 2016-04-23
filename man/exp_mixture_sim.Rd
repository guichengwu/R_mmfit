\name{exp_mixture_sim}
\alias{exp_mixture_sim}
\title{
Exponential Mixture Simulation Data
}
\description{
Thie function generates data of a mixture of two exponentially distributed random variables to supply to the gmm() function for computation and parameter estimation.
\usage{
beta_density(n, lambda1, lambda2, r)
}
\arguments{
  \item{n}{
n is the desired sample size of the simulation data.  For mixtures, we must use a very large sample size. ~> 100000
}
\item{lambda1}{Lambda1 refers to the rate of the first exponentially distributed RV.
}
\item{lambda2}{Lambda2 refers to the rate of the first exponentially distributed RV.
}  
\item{r}{r refers to the probability that the RV comes from system 1, and 1-r the probability that it comes from system 2.
}


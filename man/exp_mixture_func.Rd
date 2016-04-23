\name{exp_mixture_func}
\alias{exp_mixture_func}
\title{
 Exponential Mixture Moment Function
}
\description{
This function is supplied to the gmm() function for parameter computation when dealing with a mixture of exponentially distributed RVs.
}
\usage{
exp_mixture_func(theta, x)
}
\arguments{
  \item{theta}{
theta is a vector of the initial guesses of parameters: lambda1, lambda2, and r.    
}
  \item{x}{
x is the vector of simulation data of a mixture of two expnontially distributed random variables.
}
}

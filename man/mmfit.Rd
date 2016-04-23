\name{mmfit}
\alias{mmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mmfit function
}
\description{
The mmfit functiion will return a "mmfit" class object, which includes the estimated parameters, the standard errors and two graphs.
}
\usage{
mmfit(g, x, start, distributionType, user_density_func = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
     g is the moment function
}
  \item{x}{
     x is the simulated data
}
  \item{start}{
     start is a vector of the guess values
}
  \item{distributionType}{
     distributionType is a string that describes the distribution type. It can be "gamma" and e.t.c.
}
  \item{user_density_func}{
     User_density_func is the distribution's density function
}
}

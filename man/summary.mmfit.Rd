\name{summary.mmfit}
\alias{summary.mmfit}
\title{
   summary.mmfit function
}
\description{
     This function will summary everything the mmfit method returned, including the estimated parameters, the standard errors, and the graphs.
}
\usage{
summary.mmfit(g, x, start, distributionType, user_density_func)
}
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


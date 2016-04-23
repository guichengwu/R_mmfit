\name{plot.mmfit}
\alias{plot.mmfit}
\title{
plot.mmfit function
}
\description{
This function will draw two graphs, one is to check whether the simulated data fits the certain distribution, the other is to draw empirical cdf with its K-S confidence band. 
}
\usage{
plot.mmfit(g, x, start, distributionType, user_density_func)
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

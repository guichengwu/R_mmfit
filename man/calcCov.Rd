\name{calcCov}
\alias{calcCov}
\title{
  Analytical Covariance Calculation Function
}
\description{
This function takes the PI vector and ALPHA vector as parameters and computes the covariance of Xn and Xn+1 via mathematical covariance equations and conditional/iterative expectations.
}
\usage{
calcCov(piVector, alphaVector)
}
\arguments{
  \item{piVector}{
  The PI vector is the vector of long run probabilities that Xn is in a particular state.
}
  \item{alphaVector}{
The ALPHA vector is a uniformly random generated vector of probabilities that In will is equal to 1, thus jumping to the next higher state.
}
}



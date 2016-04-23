\name{exponential}
\alias{exponential}
\title{
Exponential RV Generation
}
\description{
 This function generates 1 exponentially distributed random variable with a probability r that it will be with rate1 or rate2.
}
\usage{
exponential(lambda1, lambda2, r)
}
\arguments{
  \item{lambda1}{
 lambda1 is the rate of the first exponential distribution to generate.
}
  \item{lambda2}{
lambda2 is the rate of the first exponential distribution to generate.
 }
  \item{r}{
r is the probability that a value will be generated with rate lambda1 and 1-r the the probabiliyt that value will be generated with rate lambda2.
}
}

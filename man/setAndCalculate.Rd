% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_utils.R
\name{setAndCalculate}
\alias{setAndCalculate}
\alias{setAndCalculateDiff}
\title{Creates a nimbleFunction for setting the values of one or more model nodes,
calculating the associated deterministic dependents and logProb values,
and returning the total sum log-probability.}
\usage{
setAndCalculate(model, targetNodes)

setAndCalculateDiff(model, targetNodes)
}
\arguments{
\item{model}{An uncompiled or compiled NIMBLE model.  This argument is required.}

\item{targetNodes}{A character vector containing the names of one or more nodes or variables in the model.  This argument is required.}
}
\description{
This nimbleFunction generator must be specialized to any model object and one or more model nodes.
A specialized instance of this nimbleFunction will set the values of the target nodes in the specified model,
calculate the associated logProbs, calculate the values of any deterministic dependents,
calculate the logProbs of any stochastic dependents,
and return the sum log-probability associated with the target nodes and all stochastic dependent nodes.
}
\details{
Calling \code{setAndCalculate(model, targetNodes)} or \code{setAndCalculate(model, targetNodes)} will return a nimbleFunction object whose \code{run} function takes a single, required argument:

targetValues: A vector of numeric values which will be put into the target nodes in the specified model object.  The length of this numeric vector much exactly match the number of target nodes.

The difference between \code{setAndCalculate} and \code{setAndCalculateDiff} is the return value of their \code{run} functions.  In the former, \code{run} returns the sum of the log probabilities of the \code{targetNodes} with the provided \code{targetValues}, while the latter returns the difference between that sum with the new \code{targetValues} and the previous values in the \code{model}.
}
\examples{
code <- nimbleCode({ for(i in 1:3) { x[i] ~ dnorm(0,1); y[i] ~ dnorm(0, 1)}})
Rmodel <- nimbleModel(code)
my_setAndCalc <- setAndCalculate(Rmodel, c('x[1]', 'x[2]', 'y[1]', 'y[2]'))
lp <- my_setAndCalc$run(c(1.2, 1.4, 7.6, 8.9))
}
\author{
Daniel Turek
}

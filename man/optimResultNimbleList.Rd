% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nimbleList_core.R
\docType{data}
\name{optimResultNimbleList}
\alias{optimResultNimbleList}
\title{Data type for the return value of \code{\link{nimOptim}}}
\format{
An object of class \code{list} of length 1.
}
\usage{
optimResultNimbleList
}
\description{
\code{\link{nimbleList}} definition for the type of \code{\link{nimbleList}} returned by \code{\link{nimOptim}}.
}
\section{Fields}{

\describe{
\item{\code{par}}{The best set of parameters found.}

\item{\code{value}}{The value of fn corresponding to par.}

\item{\code{counts}}{A two-element integer vector giving the number of calls to fn and gr respectively.}

\item{\code{convergence}}{An integer code. 0 indicates successful completion. Possible error codes are
1 indicates that the iteration limit maxit had been reached.
10 indicates degeneracy of the Nelder-Mead simplex.
51 indicates a warning from the "L-BFGS-B" method; see component message for further details.
52 indicates an error from the "L-BFGS-B" method; see component message for further details.}

\item{\code{message}}{A character string giving any additional information returned by the optimizer, or NULL.}

\item{\code{hessian}}{Only if argument hessian is true. A symmetric matrix giving an estimate of the Hessian at the solution found.}
}}

\seealso{
\code{\link{optim}}, \code{\link{nimOptim}}
}
\keyword{datasets}

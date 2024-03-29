% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cppInterfaces_nimbleFunctions.R
\name{valueInCompiledNimbleFunction}
\alias{valueInCompiledNimbleFunction}
\title{get or set value of member data from a compiled nimbleFunction using a multi-interface}
\usage{
valueInCompiledNimbleFunction(cnf, name, value)
}
\arguments{
\item{cnf}{Compiled nimbleFunction object}

\item{name}{Name of the member data}

\item{value}{If provided, the value to assign to the member data.  If omitted, the value of the member data is returned.}
}
\description{
Most nimbleFunctions written for direct user interaction allow standard R-object-like access to member data using \code{$} or \code{`[[`}.  However, sometimes compiled nimbleFunctions contained within other compiled nimbleFunctions are interfaced with a light-weight system called a multi-interface.  \code{valueInCompiledNimbleFunction} provides a way to get or set values in such cases.
}
\details{
The member data of a nimbleFunction are the objects created in \code{setup} code that are used in \code{run} code or other member functions.

Whether multi-interfaces are used for nested nimbleFunctions is controlled by the \code{buildInterfacesForCompiledNestedNimbleFunctions} option in \code{\link{nimbleOptions}}.

To see an example of a multi-interface, see \code{samplerFunctions} in a compiled MCMC interface object.
}
\author{
Perry de Valpine
}

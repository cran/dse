\name{TSmodel}
\alias{TSmodel}
\alias{TSmodel.TSmodel}
\alias{TSmodel.TSestModel}
\alias{is.TSmodel}

\title{Time Series Models}
\description{Construct or extract a "TSmodel" from objects.}
\usage{
    TSmodel(model, ...)
    \method{TSmodel}{TSmodel}(model)
    \method{TSmodel}{TSestModel}(model)
    is.TSmodel(obj)
}
\arguments{
    \item{obj}{
    An object containing an object of class TSmodel or a list containing the 
    information necessary to build an object of class "TSmodel".}
}
\details{
This is a generic method which will extract a "TSmodel" from an object (e.g. a
TSestModel. The default method will try to build an ARMA or state-space 
"TSmodel" from a list, which must contain the necessary information.

This class of objects is returned by estimation methods or can be built 
according to the description for specific sub-classes (eg "ARMA", "SS").

The \code{TSmodel} class of objects has methods for the generic functions
\code{print, test.equal, series.names, input.series.names,
output.series.names, l, roots, stability,  
forecast, feather.forecasts, horizon.forecasts,
simulate, monte.carlo.simulations}

Also, the function \code{is.TSmodel} and the functions \code{to.SS, to.ARMA, to.troll} are supported.
Other model classes inherit from the class \code{TSmodel}.

This class of objects contains a time series model. It is the class of
objects expected by many of the functions in this library.  

Sub-class (e.g. ARMA and SS for linear,
time-invariant ARMA and state space models.) are documented individually.
Many of the functions in this library are
designed for estimating and converting among various representations of
these types of models.
}

\seealso{
    \code{\link{ARMA}},
    \code{\link{SS}},
    \code{\link{TSestModel}},
    \code{\link{TSdata}}
}
%\keyword{DSE}
\keyword{ts}


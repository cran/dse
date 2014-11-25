\name{summary.estimated.models}
\alias{summary.estimated.models}
\alias{summary.TSestModel.ee}
\alias{summary.TSmodel.ee}
\alias{summary.EstEval}
\alias{summary.forecastCov}
\alias{summary.forecastCov.estimatorsWRTdata.subsets}
\alias{summary.forecastCov.estimatorsWRTtrue}
\alias{summary.MonteCarloSimulations}
\alias{summary.coef.ee}
\alias{summary.roots.ee}

\alias{print.summary.estimated.models}
\alias{print.summary.TSestModel.ee}
\alias{print.summary.TSmodel.ee}
\alias{print.summary.EstEval}
\alias{print.summary.forecastCov}
\alias{print.summary.forecastCov.estimatorsWRTdata.subsets}
\alias{print.summary.forecastCov.estimatorsWRTtrue}
\alias{print.summary.MonteCarloSimulations}
\alias{print.summary.coef.ee}
\alias{print.summary.roots.ee}

\title{Summary Specific Methods}
\description{See the generic function description.}
\usage{
    \method{summary}{estimated.models}(object)
    \method{summary}{TSestModel.ee}(object)
    \method{summary}{TSmodel.ee}(object)
    \method{summary}{EstEval}(object)
    \method{summary}{forecastCov}(object, horizons, series)
    \method{summary}{forecastCov.estimatorsWRTdata.subsets}(object)
    \method{summary}{forecastCov.estimatorsWRTtrue}(object, digits)
    \method{summary}{MonteCarloSimulations}(object, series, periods)
    \method{summary}{coef.ee}(object, verbose=TRUE)
    \method{summary}{roots.ee}(object, verbose=TRUE)
    
    \method{print}{summary.estimated.models}(x, digits)
    \method{print}{summary.TSestModel.ee}(x, digits)
    \method{print}{summary.TSmodel.ee}(x, digits)
    \method{print}{summary.EstEval}(x, digits)
    \method{print}{summary.forecastCov}(x, digits)
    \method{print}{summary.forecastCov.estimatorsWRTdata.subsets}(x, digits)
    \method{print}{summary.forecastCov.estimatorsWRTtrue}(x, digits)
    \method{print}{summary.MonteCarloSimulations}(x, digits)
    \method{print}{summary.coef.ee}(x, digits)
    \method{print}{summary.roots.ee}(x, digits)
}
\arguments{
    \item{object}{an object for which a summary is to be printed.}
    \item{x}{an object for which a summary is to be printed.}
    \item{digits}{a non-null value is used to indicate the number of significant
       digits. If \code{digits} is \code{NULL} then the value of digits 
       specified by \code{options} is used.}
    \item{horizons}{optional integer vector indicating horizons at which the 
       summary should be calculated.}
    \item{series}{
       The series which should be plotted. The default NULL gives all series.}
    \item{periods}{optional integer vector indicating periods at which the 
       summary should be calculated.}
    \item{verbose}{logical indicating if a longer summary should be produced.}
}

\seealso{
    \code{\link{summary}}
    \code{\link{print}}
}
%\keyword{DSE}
\keyword{ts}



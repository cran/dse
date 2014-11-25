\name{summary.estimated.models}
\alias{summary.estimated.models}
\alias{summary.TSestModel.ee}
\alias{summary.TSmodel.ee}
\alias{summary.EstEval}
\alias{summary.forecastCov}
\alias{summary.forecastCovEstimatorsWRTdata.subsets}
\alias{summary.forecastCovEstimatorsWRTtrue}
\alias{summary.MonteCarloSimulations}
\alias{summary.coef.ee}
\alias{summary.roots.ee}

\alias{print.summary.estimated.models}
\alias{print.summary.TSestModel.ee}
\alias{print.summary.TSmodel.ee}
\alias{print.summary.EstEval}
\alias{print.summary.forecastCov}
\alias{print.summary.forecastCovEstimatorsWRTdata.subsets}
\alias{print.summary.forecastCovEstimatorsWRTtrue}
\alias{print.summary.MonteCarloSimulations}
\alias{print.summary.coef.ee}
\alias{print.summary.roots.ee}

\title{Summary Specific Methods}
\description{See the generic function description.}
\usage{
    \method{summary}{estimated.models}(object, ...)
    \method{summary}{TSestModel.ee}(object, ...)
    \method{summary}{TSmodel.ee}(object, ...)
    \method{summary}{EstEval}(object, ...)
    \method{summary}{forecastCov}(object, horizons=object$horizons,
        series=seq(nseriesOutput(object$data)), ...)
    \method{summary}{forecastCovEstimatorsWRTdata.subsets}(object, ...)
    \method{summary}{forecastCovEstimatorsWRTtrue}(object,
        digits=options()$digits, ...)
    \method{summary}{MonteCarloSimulations}(object, series=NULL, periods=1:3, ...)
    \method{summary}{coef.ee}(object, verbose=TRUE, ...)
    \method{summary}{roots.ee}(object, verbose=TRUE, ...)
    
    \method{print}{summary.estimated.models}(x, digits=options()$digits, ...)
    \method{print}{summary.TSestModel.ee}(x, digits=options()$digits, ...)
    \method{print}{summary.TSmodel.ee}(x, digits=options()$digits, ...)
    \method{print}{summary.EstEval}(x, digits=options()$digits, ...)
    \method{print}{summary.forecastCov}(x, digits=options()$digits, ...)
    \method{print}{summary.forecastCovEstimatorsWRTdata.subsets}(x,
       digits=options()$digits, ...)
    \method{print}{summary.forecastCovEstimatorsWRTtrue}(x,
       digits=options()$digits, ...)
    \method{print}{summary.MonteCarloSimulations}(x, digits=options()$digits, ...)
    \method{print}{summary.coef.ee}(x, digits=options()$digits, ...)
    \method{print}{summary.roots.ee}(x, digits=options()$digits, ...)
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
    \item{...}{arguments passed to other methods.}
}

\seealso{
    \code{\link{summary}}
    \code{\link{print}}
}
%\keyword{DSE}
\keyword{ts}



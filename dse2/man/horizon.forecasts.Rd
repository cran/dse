\name{horizon.forecasts}
\alias{horizon.forecasts}
\alias{horizon.forecasts.TSdata}
\alias{horizon.forecasts.TSmodel}
\alias{horizon.forecasts.TSestModel}
\alias{horizon.forecasts.compiled}
\alias{horizon.forecasts.compiled.ARMA}
\alias{horizon.forecasts.compiled.SS}
\alias{horizon.forecasts.forecast.cov}
\alias{is.horizon.forecasts}

\title{Calculate forecasts at specified horizons}
\description{Calculate forecasts at specified horizons.}
\usage{
	horizon.forecasts(model, data, horizons=1:4, 
	    discard.before=minimum.startup.lag(model), compiled=.DSECOMPILED)
	is.horizon.forecasts(obj)
}
\arguments{
\item{ model}{An object of class TSmodel.}
\item{ data}{An object of class TSdata} 
\item{ horizons}
{A vector of integers indicating the horizon at which forecasts should be 
produced.}
\item{ discard.before}{Period before which forecasts are not calculated.}
\item{ compiled}{If T compiled code is called.}
}
\value{
The result is a list of class horizon.forecasts with elements model (a 
TSmodel), data, horizons, discard.before, and horizon.forecasts.
horizon.forecasts is an array with three dimension: 
c(length(horizons), dim(model$data)).
Projections are not calculated before discard.before or after
the end of output.data(data).
Each horizon is aligned so that horizon.forecasts[h,t,] contains the forecast
for the data point output.data(data)[t,] (from horizon[h] periods prior).
}
\details{
Calculate multiple 'horizon'-step ahead forecasts 
ie. calculate forecasts but return only those indicated by horizons.
   Thus, for example, the result of
 horizon.forecasts(model, data horizons=c(1,5))    
 would be the one-step ahead and five step ahead forecasts.
}
\seealso{
\code{\link{feather.forecasts}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <-  horizon.forecasts(model, eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


\name{feather.forecasts}
\alias{feather.forecasts}
\alias{feather.forecasts.TSdata}
\alias{feather.forecasts.TSmodel}
\alias{feather.forecasts.TSestModel}
\alias{is.feather.forecasts}

\title{Multiple Horizon-Step Ahead Forecasts}
\description{Calculate multiple horizon-step ahead forecasts.}
\usage{
	feather.forecasts(model, data, 
	   from.periods = 10*seq(periods(data)/10), horizon=36, ...)
	is.feather.forecasts(obj)
}
\arguments{
\item{ model}{An object of class TSmodel.}
\item{ data}{An object of class TSdata.}
\item{ from.periods.}{The starting points to use for forecasts}
\item{ horizon}{The number of periods to forecast.}
\item{ ...}{For a TSmodel additional arguments are passed to l()}
}
\value{
The result is a list of class feather.forecasts with elements $model (a 
TSestModel), $data, $from.periods, $feather.forecasts.
$feather.forecasts is a list with length(from.periods) elements, each of which
is a tframed matrix.
There is a plot method for this class.
}
\details{
Calculate multiple horizon-step ahead forecasts 
ie. use the samples indicated by from.periods to calculate forecasts for 
horizon periods.
Thus, for example, the result of
feather.forecasts(model, data, from.periods=c(200,250,300))
would be forecasts for 1 through 36 steps ahead (the default),
starting at the 200th,250th, and 300th point of output.data(data).
This function assumes that input.data(data) (the exogenous 
variable) is as long as necessary for the most future forecast.
}
\seealso{
\code{\link{forecast}},
\code{\link{horizon.forecasts}}
}
\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
model <- est.VARX.ls(egJofF.1dec93.data)
pr <- feather.forecasts(model, egJofF.1dec93.data)
tfplot(pr)
}
\keyword{DSE}
\keyword{ts}


\name{forecast}
\alias{forecast}
\alias{forecast.TSdata}
\alias{forecast.TSmodel}
\alias{forecast.TSestModel}
\alias{is.forecast}

\title{Forecast Multiple Steps Ahead}
\description{Calculate forecasts multiple steps ahead.}
\usage{
	forecast(model, data,  horizon=36, percent=100, 
		conditioning.inputs=NULL, conditioning.inputs.forecasts=NULL)
	is.forecast.cov(x)
}
\arguments{
\item{ model}{An object of class TSmodel.}
\item{ data}{An object of class TSdata.}
\item{ conditioning.inputs}
{A time series matrix or list of time series matrices to use as input variables.}
\item{ conditioning.inputs.forecasts}
{A time series matrix or list of time series matrices to append to input
variables for the forecast periods.}
\item{ horizon }
{The number of periods to forecast.}
\item{ percent}
{A vector indication percentages of the last input to use for forecast periods.}
}
\value{
The result is an object of class forecast which is a list with 
elements \code{model}, \code{horizon}, \code{conditioning.inputs}, 
\code{percent}, \code{pred} and \code{forecast}.
The element \code{forecast} is a list with TSdata objects as elements, one for each element
in the list conditioning.inputs.
The element \code{pred} contains the one-step ahead forecasts for the preiods when output data
is available. 
There is a plot method for this class.
}
\details{
Calculate (multiple) forecasts from the end of data to a horizon determined
either from supplied input data or the argument horizon (more details below).
In  the case of a model with no inputs the horizon is determined by the
argument horizon.
In the case of models with inputs, on which the forecasts
are conditioned, the argument horizon is ignored (except when percent is
specified) and the actual horizon is determined by the inputs in the 
following way:
If inputs are not specified by optional arguments (as below) then the default
will be to use input.data(data). This will be the same as the function l() unless
input.data(data) is longer than
output.data(data) (after NAs are trimmed from each separately).
Otherwise, if conditioning.inputs is specified it is used for input.data(data).
It must be a time series matrix or a list of time series matrices each
of which is used in turn as input.data(data). The default above is the same as
forecast(model, trim.na(data), conditioning.inputs=trim.na(input.data(data)) )
Otherwise, if conditioning.inputs.forecasts is specified it is appended 
to input.data(data). It must be a time series  
matrix or a list of time series matrices each of which is 
appended to input.data(data) and the concatenation used as conditioning.inputs.
Both conditioning.inputs and conditioning.inputs.forecasts should not be
specified.
Otherwise, if percent is specified then conditioning.inputs.forecasts are set
to percent/100 times the value of input corresponding to the last period
of output.data(data) and used for horizon periods. percent can be a vector, 
in which case each value is applied in turn. ie c(90,100,110) would would 
give results for conditioning.input.forecasts 10 percent above and below 
the last value of input.
}
\seealso{
\code{\link{feather.forecasts}},
\code{\link{horizon.forecasts}}
}
\examples{
	if(is.R()) data("egJofF.1dec93.data", package="dse1")
	model <- est.VARX.ls(window(egJofF.1dec93.data, end=c(1985,12)))
	pr <- forecast(model, conditioning.inputs=input.data(egJofF.1dec93.data))
	#tfplot(pr) Rbug 0.90.1
	is.forecast(pr)
} 
\keyword{DSE}
\keyword{ts}


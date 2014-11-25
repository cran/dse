\name{forecasts}
\alias{forecasts}
\alias{forecasts.forecast}
\alias{forecasts.feather.forecasts}
\alias{forecasts.horizon.forecasts}
\title{Extract Forecasts}
\description{Extract forecasts from and object.}
\usage{
	forecasts(obj, ...)
}
\arguments {
	\item{obj}{An object which contains forecasts.}
}
\value{The forecasts from an object which contains forecasts.}
\details{
	This generic method extracts the forecasts (only) from objects returned 
	by other methods that calculate forecasts. Usually the objects returned
	by the methods which calculate forecasts contain additional
	information which is not returned by this extractor. 
}
\seealso{
	\code{\link{forecast}}
}
\examples{
	if(is.R()) data("egJofF.1dec93.data", package="dse1")
	model <- est.VARX.ls(window(egJofF.1dec93.data, end=c(1985,12)))
	pr <- forecast(model, conditioning.inputs=input.data(egJofF.1dec93.data))
	z <- forecasts(pr)
} 
\keyword{DSE}
\keyword{ts}


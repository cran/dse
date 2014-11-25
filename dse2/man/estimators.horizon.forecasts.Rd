\name{estimators.horizon.forecasts}
\alias{estimators.horizon.forecasts}
\title{Estimate models and forecast at given horizons}
\description{Estimate models and forecast at given horizons.}
\usage{
	estimators.horizon.forecasts(data, 
                       estimation.sample=.5, horizons=1:12,quiet=F,
                       estimation.methods=NULL)
}
\arguments{
	\item{data}{A TSdata object.}
	\item{estimation.methods}{
	A list of estimation methods to use. (See estimate.models.)}
	\item{estimation.sample}{The portion  of the sample to use for estimation.}
	\item{horizons}{The horizons for which forecasts are to be produced.}
	\item{quiet}{If true no estimation information is printed.}
}
\value{A list of forecasts at different horizons as returned by horizon.forecasts.
}
\details{
estimation.sample indicates the part of the data to use for estimation.
If estimation.sample is less than or equal 1.0 it is
used to indicate the portion of points to use for estimation.
Otherwise it should be an integer and is used to indicate the number
of points from the beginning of the sample to use for estimation. 
}
\seealso{
\code{\link{estimate.models}},
\code{\link{horizon.forecasts}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  estimators.horizon.forecasts(eg1.DSE.data.diff, 
	estimation.methods=list(est.VARX.ls=list(max.lag=3), 
				est.VARX.ar=list(max.lag=3)))
}
\keyword{DSE}
\keyword{ts}


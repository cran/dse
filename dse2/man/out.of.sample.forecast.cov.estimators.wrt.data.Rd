\name{out.of.sample.forecast.cov.estimators.wrt.data}
\alias{out.of.sample.forecast.cov.estimators.wrt.data}
\title{Calculate Out-of-Sample Forecasts}
\description{Calculate out-of-sample forecasts.}
\usage{
	out.of.sample.forecast.cov.estimators.wrt.data(data, zero=F, trend=F,
                       estimation.sample=.5, horizons=1:12,quiet=F,
                       estimation.methods=NULL)

}
\arguments{
\item{ data}{An object of class TSdata.}
\item{ estimation.methods}{A list as used by estimate.models.}
\item{ zero}{If T then forecast.cov is also calculated for a forecast of zero.}
\item{ trend}
{If T then forecast.cov is also calculated for a forecast of a linear trend.}
\item{ estimation.sample}
{Indicates the portion of the data to use for estimation.
If estimation.sample is an integer then it is used to indicate the number
of points in the sample to use for estimation. If it is a fracton it is
used to indicate the portion of points to use for estimation. The remainder
of the sample is used for evaluating forecasts.}
\item{ horizons}{Horizons for which forecast covariance should be calculated.}
\item{ quiet}{If T then estimation information is not printed.}
}
\value{An object as returned by forecast.cov.estimators.wrt.data.}
\details{
The data is slpit into a sub-sample used for estimation and another
sub-sample used for calculating the forecast covariance.
}
\seealso{
\code{\link{forecast.cov.estimators.wrt.data}},
\code{\link{estimate.models}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  out.of.sample.forecast.cov.estimators.wrt.data(eg1.DSE.data.diff,
           estimation.methods=list(est.VARX.ls=list(max.lag=4)))
}
\keyword{DSE}
\keyword{ts}


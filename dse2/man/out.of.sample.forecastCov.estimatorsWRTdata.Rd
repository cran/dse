\name{out.of.sample.forecastCov.estimatorsWRTdata}
\alias{out.of.sample.forecastCov.estimatorsWRTdata}
\title{Calculate Out-of-Sample Forecasts}
\description{Calculate out-of-sample forecasts.}
\usage{
    out.of.sample.forecastCov.estimatorsWRTdata(data, zero=F, trend=F,
                       estimation.sample=.5, horizons=1:12,quiet=F,
                       estimation.methods=NULL, compiled=.DSECOMPILED)
}
\arguments{
    \item{data}{an object of class TSdata.}
    \item{estimation.methods}{a list as used by estimate.models.}
    \item{zero}{if T then forecastCov is also calculated for a forecast of zero.}
    \item{trend}{
    if T then forecastCov is also calculated for a forecast of a linear trend.}
    \item{estimation.sample}{
    indicates the portion of the data to use for estimation.
If estimation.sample is an integer then it is used to indicate the number
of points in the sample to use for estimation. If it is a fracton it is
used to indicate the portion of points to use for estimation. The remainder
of the sample is used for evaluating forecasts.}
    \item{horizons}{horizons for which forecast covariance should be calculated.}
    \item{quiet}{if T then estimation information is not printed.}
    \item{compiled}{a logical indicating if compiled code should be
      used. (Usually true except for debugging.)}
}
\value{An object as returned by forecastCov.estimatorsWRTdata.}
\details{
The data is slpit into a sub-sample used for estimation and another
sub-sample used for calculating the forecast covariance.
}
\seealso{
\code{\link{forecastCov.estimatorsWRTdata}},
\code{\link{forecastCov.estimatorsWRTtrue}},
\code{\link{estimate.models}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  out.of.sample.forecastCov.estimatorsWRTdata(eg1.DSE.data.diff,
           estimation.methods=list(est.VARX.ls=list(max.lag=4)))
}
%\keyword{DSE}
\keyword{ts}


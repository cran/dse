\name{forecast.cov.estimators.wrt.data}
\alias{forecast.cov.estimators.wrt.data}
\alias{is.forecast.cov.estimators.wrt.data}

\title{Calculate Forecast Cov of Estimators WRT Data}
\description{forecast covariance of estimated models wrt given sample}
\usage{
    forecast.cov.estimators.wrt.data(data, estimation.sample=NULL, 
                       compiled=.DSECOMPILED, discard.before=10,
                       horizons=1:12, zero=F, trend=F,quiet=F,
                       estimation.methods=NULL)
    is.forecast.cov.estimators.wrt.data(obj)
}
\arguments{
    \item{data}{An object of class TSdata.}
    \item{estimation.methods}{A list as used by estimate.models.}
    \item{discard.before}{
    An integer indicating the number of points in the
    beginning of forecasts to discard for calculating covariances.}
    \item{zero}{If T then forecast.cov is also calculated for a forecast of zero.}
    \item{trend}{
    If T then forecast.cov is also calculated for a forecast of a linear trend.}
    \item{estimation.sample}{
    An integer indicating the number of points in the
sample to use for estimation. If it is NULL the whole sample is used.}
    \item{horizons}{Horizons for which forecast covariance should be calculated.}
    \item{quiet}{If T then estimation information is not printed.}
    \item{horizons}{
    Horizons for which forecast covariance should be calculated.}
}
\value{
A list with the forecast covariance for supplied models on the given
sample. This is in the element \code{forecast.cov} of the result. Other
elements contain information in the arguments. 
}
\details{
Calculate the forecasts cov of models estimated from data with estimation
methods indicated by estimation.methods  (see estimate.models).
estimation.sample is an integer indicating the number of points in the
sample to use for estimation. If it is NULL the whole sample is used.
}
\seealso{
\code{\link{out.of.sample.forecast.cov.estimators.wrt.data}},
\code{\link{estimate.models}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- forecast.cov.estimators.wrt.data(eg1.DSE.data.diff, 
    estimation.methods=list(est.VARX.ls=list(max.lag=4)))
}
%\keyword{DSE}
\keyword{ts}


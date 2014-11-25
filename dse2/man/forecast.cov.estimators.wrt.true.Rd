\name{forecast.cov.estimators.wrt.true}
\alias{forecast.cov.estimators.wrt.true}
\alias{is.forecast.cov.estimators.wrt.true}

\title{Compare Forecasts Cov Relative to True Model Output}
\description{}
\usage{
	forecast.cov.estimators.wrt.true(true.model, Spawn=.SPAWN, rng=NULL,
                       simulation.args=NULL,
                       est.replications=2, pred.replications=2,
                       discard.before=10, horizons=1:12,quiet=F,
                       estimation.methods=NULL)
	is.forecast.cov.estimators.wrt.true(obj)
}
\arguments{
	\item{true.model}{An object of class TSmodel.}
	\item{estimation.methods}{A list as used by estimate.models.}
	\item{simulation.args}{An arguments to be passed to simulate.}
	\item{discard.before}{
	An integer indicating the number of points in the
	beginning of forecasts to discard for calculating covariances.}
	\item{zero}{If T then forecast.cov is also calculated for a forecast of zero.}
	\item{trend}{
	If T then forecast.cov is also calculated for a forecast of a linear trend.}
	\item{horizons}{Horizons for which forecast covariance should be calculated.}
	\item{horizons}{Horizons for which forecast covariance should be calculated.}
	\item{rng}{If specified then it is used to set RNG.}
	\item{Spawn}{If T then Splus For loops are used.}
	\item{quiet}{If T then some messages are not printed.}
}
\value{
The returned results has element
  \code{forecast.cov.true, forecast.cov.zero, forecast.cov.trend} containing 
     covariances averaged over estimation replications and simulation
     replications (forecasts will not change but simulated data will).
   \code{forecast.cov} a list of the same length as estimation.methods with each
     element containing covariances averaged over estimation replications 
     and simulation replications.
   \code{estimated.models} a list of length est.replications, with each elements as
     returned by estimate.models, thus each element has \code{multi.model} as a
     subelement containing models for different estimation techniques.  
     So, eg. \code{estimated.models[[2]]$multi.model[[1]]}  in the result will
     be the model from the first estimation technique in the second replication.
}
\details{
Calculate the forecasts cov of models estimated from simulations of 
  true.model with estimation methods indicated by estimation.methods (see 
       estimate.models). 
This function makes multiple calls to forecast.cov.wrt.true.
}
\seealso{
\code{\link{forecast.cov.wrt.true}}
\code{\link{forecast.cov.estimators.wrt.data}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
true.model <- est.VARX.ls(eg1.DSE.data.diff) # just to have a starting model
z <-  forecast.cov.estimators.wrt.true(true.model, 
	estimation.methods=list(est.VARX.ls=list(max.lag=4)))
}
\keyword{DSE}
\keyword{ts}


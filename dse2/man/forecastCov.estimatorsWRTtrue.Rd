\name{forecastCov.estimatorsWRTtrue}
\alias{forecastCov.estimatorsWRTtrue}
\alias{is.forecastCov.estimatorsWRTtrue}

\title{Compare Forecasts Cov Relative to True Model Output}
\description{}
\usage{
    forecastCov.estimatorsWRTtrue(true.model, Spawn=.SPAWN, rng=NULL,
                       simulation.args=NULL,
                       est.replications = 2, pred.replications = 2,
                       discard.before = 10, horizons = 1:12, quiet =FALSE,
                       estimation.methods=NULL, compiled=.DSECOMPILED)
    is.forecastCov.estimatorsWRTtrue(obj)
}
\arguments{
    \item{true.model}{An object of class TSmodel.}
    \item{estimation.methods}{A list as used by estimate.models.}
    \item{simulation.args}{an arguments to be passed to simulate.}
    \item{est.replications}{An arguments to be passed to simulate.}
    \item{pred.replications}{An arguments to be passed to simulate.}
    \item{discard.before}{
      An integer indicating the number of points in the
      beginning of forecasts to discard for calculating covariances.}
    \item{zero}{If TRUE then forecastCov is also calculated for a forecast of zero.}
    \item{trend}{If TRUE then forecastCov is also calculated for a 
      forecast of a linear trend.}
    \item{horizons}{Horizons for which forecast covariance should be calculated.}
    \item{rng}{If specified then it is used to set RNG.}
    \item{Spawn}{If TRUE then Splus For loops are used.}
    \item{quiet}{If TRUE then some messages are not printed.}
    \item{compiled}{a logical indicating if the compiled version of the code
       should be used. (FALSE would typically only be used for debugging.)}
    \item{obj}{an object.}
}
\value{
The returned results has element
  \code{forecastCov.true, forecastCov.zero, forecastCov.trend} containing 
     covariances averaged over estimation replications and simulation
     replications (forecasts will not change but simulated data will).
   \code{forecastCov} a list of the same length as estimation.methods with each
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
This function makes multiple calls to forecastCovWRTtrue.
}
\seealso{
\code{\link{forecastCovWRTtrue}}
\code{\link{forecastCov.estimatorsWRTdata}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
true.model <- est.VARX.ls(eg1.DSE.data.diff) # just to have a starting model
z <-  forecastCov.estimatorsWRTtrue(true.model, 
    estimation.methods=list(est.VARX.ls=list(max.lag=4)))
}
%\keyword{DSE}
\keyword{ts}


\name{select.forecast.cov}
\alias{select.forecast.cov}
\title{Select Forecast Covariances Meeting Criteria}
\description{Select forecast covariances meeting given criteria.}
\usage{
    select.forecast.cov(obj, select.series=1, 
    select.cov.best=1,
    select.cov.bound=NULL,
    ranked.on.cov.bound=NULL,
    verbose=T)
}
\arguments{
    \item{obj}{An object as returned by mine.strip.}
    \item{select.series}{
    An indication of series to which the tests should be applied.}
    \item{select.cov.best}{The number of 'best' forecasts to select.}
    \item{select.cov.bound}{A bound to use as criteria for selection.}
    \item{verbose}{If verbose=T then summary results are printed.}
}
\value{
The returned result is a forecast.cov object like obj, but filtered
to remove models which do not meet criteria.
}
\details{
Select models with forecast covariance for select.series meeting criteria.
The default select.cov.best=1 selects the best model at each horizon.
select.cov.best=3 would select the best 3 models at each horizon.
If select.cov.bound is not NULL then  select.cov.best is ignored and
any model which is better than the bound at all horizons is selected.
select.cov.bound can be a vector of the same length as select.series,
in which case corresponding elements are applied to the different series.
}
\seealso{
\code{\link{min.forecast.cov}},
\code{\link{exclude.forecast.cov}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
z <-  select.forecast.cov(z)
tfplot(select.forecast.cov(z, select.cov.bound=20000))
tfplot(select.forecast.cov(z, select.cov.best=1))
}
%\keyword{DSE}
\keyword{ts}


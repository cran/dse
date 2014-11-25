\name{exclude.forecast.cov}
\alias{exclude.forecast.cov}
\title{Filter Object to Remove Forecasts}
\description{Filter object to remove forecasts.}
\usage{
    exclude.forecast.cov(obj, exclude.series=NULL)
}
\arguments{
    \item{obj}{An object as returned by mine.strip.}
    \item{exclude.series}{An indication of series to which should be excluded.}
}
\value{
The returned result is a forecast.cov object like obj, but filtered to
remove any forecasts from models which depend on the series which are
indicated for exclusion.
}
\details{
Exclude results which depend on the indicated series from a 
(forecast.cov.estimators.wrt.data.subsets forecast.cov) object.
}
\seealso{
\code{\link{min.forecast.cov}},
\code{\link{select.forecast.cov}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
z <-  exclude.forecast.cov(z, exclude.series=3)
}
%\keyword{DSE}
\keyword{ts}


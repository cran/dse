\name{exclude.forecastCov}
\alias{exclude.forecastCov}
\title{Filter Object to Remove Forecasts}
\description{Filter object to remove forecasts.}
\usage{
    exclude.forecastCov(obj, exclude.series=NULL)
}
\arguments{
    \item{obj}{An object as returned by mine.strip.}
    \item{exclude.series}{An indication of series to which should be excluded.}
}
\value{
The returned result is a forecastCov object like obj, but filtered to
remove any forecasts from models which depend on the series which are
indicated for exclusion.
}
\details{
Exclude results which depend on the indicated series from a 
(forecastCovEstimatorsWRTdata.subsets forecastCov) object.
}
\seealso{
\code{\link{min.forecastCov}},
\code{\link{select.forecastCov}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
z <-  exclude.forecastCov(z, exclude.series=3)
}
%\keyword{DSE}
\keyword{ts}



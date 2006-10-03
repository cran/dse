\name{ytoypc}
\alias{ytoypc}
\alias{ytoypc.default}
\title{Convert to year to year percent change}
\description{Convert level data to year over year percent change.}
\usage{
    ytoypc(obj, names = paste("y to y \%ch", seriesNames(obj)))
    \method{ytoypc}{default}(obj, names = paste("y to y \%ch", seriesNames(obj)))
}
\arguments{
    \item{obj}{A time series object.}
    \item{names}{names for the new series (but is details).}
}
\value{A time series vector or matrix.}
\details{
  The result is a time series of the year over year percent change. This uses 
percentChange with lag=frequency(obj).

  The \code{names} are not applied to the new series if the global option
  ModSeriesNames is FALSE. This can be set 
  with \code{options(ModSeriesNames=FALSE)}. This provides a convenient 
  mechanism to prevent changing series labels on plot axis, when the title
  may indicate that data is in year-to-year percent change so the axis label
  does not need this.
}
\seealso{
\code{\link{percentChange}}
}
\examples{
z <- matrix(100 + rnorm(200),100,2)
z[z == 0] <- 1 # not to likely, but it can happen
z <- ytoypc(z)
}
\keyword{ts}


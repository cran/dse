\name{tframeTime}
\alias{time.tf}
\alias{time.tframe}
\alias{tframeTime}
\alias{tframeTime.default}
\alias{tframeTime.tftframe}
\alias{tframeTime.tstframe}

\title{Tframe time}
\description{Extract time from a tframe (coresponding to each data point).}
\usage{
    \method{time}{tf}(x)
    \method{time}{tframe}(x)

    tframeTime(x)
    \method{tframeTime}{default}(x)
    \method{tframeTime}{tftframe}(x)
    \method{tframeTime}{tstframe}(x)
    }

\arguments{
    \item{x}{A tframe or tframed object.}
}
\value{A time series of the dates of observations of x.}
\details{
    This returns the time (date) of each observation, as might be used on
    the axis of a graph or lables of a printout.
}
\seealso{
    \code{\link{time}}
    \code{\link{periods}}
    \code{\link{tframe}}
}
\examples{
    z <- tframed(matrix(rnorm(200), 100,2), 
            tf=list(start=c(1982,1), frequency=12))
    time(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



\name{tfdiff}

\alias{tfdiff}
\alias{tfdiff.default}
\alias{tfdiff.tframe}
\alias{diff.tframed}
%\alias{lag}
%\alias{lag.default}
\title{Time Series Differencing}
\description{Difference a tframed object.}
\usage{
    tfdiff(x, lag=1, differences=1)
    \method{tfdiff}{default}(x,lag=1, differences=1)
    \method{tfdiff}{tframe}(x,lag=1, differences=1)
    \method{diff}{tframed}(x,lag=1, differences=1, ...)
}
\arguments{
    \item{x}{a tframed object.}
    \item{lag}{difference calculated relative to lag periods previous.}
    \item{differences}{order of differencing.}
    \item{...}{arguments to be passed to other methods.}
}
\details{A time framed object is created by differencing the number of times
   indicated by differences at a lagged number of periods indicated by lag. The
   default is take the difference from data one period previous.
}
\seealso{
    \code{\link[base]{diff}},
    \code{\link[ts]{lag}}
}
\examples{
    z <- ts(rnorm(100), start=c(1982,1), frequency=12)
    tfstart(z)
    tfperiods(z)
    z <- tfdiff(z)
    tfstart(z)
    tfperiods(z)
}
\keyword{internal}



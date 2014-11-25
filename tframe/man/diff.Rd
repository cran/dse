\name{diff}

\alias{diff}
\alias{diff.tframe}
\alias{lag}
\alias{lag.default}
\title{Time Series Shifting Operations}
\description{lag a tframed object.}
\usage{
    diff(x, ...)
    \method{diff}{tframe}(x, lag = 1, differences = 1, ...)

%    lag(x, ...)
%    \method{lag}{default}(x, ...)
}
\arguments{
    \item{x}{a tframed object.}
    \item{lag}{difference calculated relative to lag periods previous.}
    \item{differences}{order of differencing.}
    \item{...}{arguments passed to other methods.}
}
\details{A time framed object is created by differencing the number of times
   indicated by differences at a lagged number of periods indicated by lag. The
   default is take the difference from data one period previous.
   For lag see documentation for the generic (from library ts in R).
}
\seealso{
    \code{\link[base]{diff}}
    \code{\link[ts]{lag}}
}
\examples{
    z <- ts(rnorm(100), start=c(1982,1), frequency=12)
    start(z)
    periods(z)
    z <- diff(z)
    start(z)
    periods(z)
}
\keyword{internal}



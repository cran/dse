\name{check.tframeConsistent}
\alias{check.tframeConsistent}
\alias{check.tframeConsistent.default}
\alias{check.tframeConsistent.stamped}

\title{Check for a Consistent tframe}
\description{
    Check if tframe and a time series are consistent with one another.}
\usage{
    check.tframeConsistent(tf, x)
    \method{check.tframeConsistent}{default}(tf, x)
    \method{check.tframeConsistent}{stamped}(tf, x)
    }

\arguments{
    \item{tf}{A tframe)}
    \item{x}{An object)}
}
\value{A logical scalar.}
\details{Check if the number of periods in the tframe corresponds to the
    number of observations in the time series.
}
\seealso{
    \code{\link{is.tframe}}
}

\examples{
    z <- ts(rnorm(100), start=c(1982,1), frequency=12)
    check.tframeConsistent(tframe(z), rnorm(100))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



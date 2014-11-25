\name{periods}

\alias{start}
\alias{start.tframe}
\alias{start.tstframe}

\alias{end}
\alias{end.tframe}
\alias{end.tstframe}

\alias{frequency}
\alias{frequency.tframe}

\alias{periods}
\alias{periods.default}
\alias{periods.tframe}
\alias{periods.stamped}

\alias{time}
\alias{time.tframe}

\alias{diff}
\alias{diff.tframe}

\title{Generic Time Series Parameters}
\description{
    Functions for extracting time information from tframes or 
    tframed objects.
}

\usage{
    start(x, ...)
    \method{start}{tframe}(x, ...)
    \method{start}{tstframe}(x, ...)

    end(x, ...)
    \method{end}{tframe}(x, ...)
    \method{end}{tstframe}(x, ...)

    frequency(x, ...)
    \method{frequency}{tframe}(x, ...)

    periods(x)
    \method{periods}{default}(x)
    \method{periods}{tframe}(x)
    \method{periods}{stamped}(x)

    time(x, ...)
    \method{time}{tframe}(x, ...)

    diff(x, ...)
    \method{diff}{tframe}(x, lag = 1, differences = 1, ...)
}

\arguments{
    \item{x}{a tframe or a tframed object.}
    \item{lag}{difference calculated relative to lag periods previous.}
    \item{differences}{order of differencing.}
    \item{...}{arguments to be passed to other methods.}
}
\value{Depends}
\details{
    The methods start and end return the start or end date of a tframe or
    tframed object. Periods return the number of observations (time points).
    frequency returns the frequency of observation, 
    typically the number of observations in a year for economic data, but
    possibly something else in other contexts.  The concept of frequency is
    not very consistently defined for time series data, and the use of 
    the frequency method should probably be avoided where possible. In
    practice it seems rarely necessary, but the method makes porting of 
    older code much easier. 
}

\seealso{
    \code{\link{tframe}}, 
    \code{\link{tframed}}
    \code{\link[base]{start}}
    \code{\link[base]{end}}
    \code{\link[base]{frequency}}
    \code{\link[base]{periods}}
    \code{\link[base]{time}}
    \code{\link[base]{diff}}
}
\examples{
    z <- ts(rnorm(100), start=c(1982,1), frequency=12)
    start(z)
    z <- tframed(matrix(rnorm(200), 100,2), 
            tf=list(start=c(1982,1), frequency=12))
    end(z)
    periods(z)
    time(z)
}

\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}


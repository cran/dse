\name{tframeStart}
\alias{start.tf}
\alias{start.tframe}
\alias{tframeStart}
\alias{tframeStart.default}
\alias{tframeStart.tftframe}
\alias{tframeStart.tstframe}
\alias{end.tf}
\alias{end.tframe}
\alias{tframeEnd}
\alias{tframeEnd.default}
\alias{tframeEnd.tftframe}
\alias{tframeEnd.tstframe}
\alias{frequency.tf}
\alias{frequency.tframe}
\alias{tframeFrequency}
\alias{tframeFrequency.default}
\alias{tframeFrequency.tftframe}
\alias{tframeFrequency.tstframe}
\alias{periods}
\alias{periods.default}
\alias{periods.tf}
\alias{periods.tframe}
\alias{tframePeriods}
\alias{tframePeriods.default}
\alias{tframePeriods.stamped}
\alias{tframePeriods.tftframe}
\alias{tframePeriods.tstframe}
\alias{periods.tsp}

\title{Generic Time Series Parameters}
\description{
    Functions for extracting time information from tframes or 
    tframed objects.
}

\usage{
    \method{start}{tf}(x)
    \method{start}{tframe}(x)
    tframeStart(x)
    \method{tframeStart}{default}(x)
    \method{tframeStart}{tftframe}(x)
    \method{tframeStart}{tstframe}(x)

    \method{end}{tf}(x)
    \method{end}{tframe}(x)
    tframeEnd(x)
    \method{tframeEnd}{default}(x)
    \method{tframeEnd}{tftframe}(x)
    \method{tframeEnd}{tstframe}(x)

    \method{frequency}{tf}(x)
    \method{frequency}{tframe}(x)
    tframeFrequency(x)
    \method{tframeFrequency}{default}(x)
    \method{tframeFrequency}{tftframe}(x)
    \method{tframeFrequency}{tstframe}(x)

    periods(x)
    \method{periods}{default}(x)
    \method{periods}{tf}(x)
    \method{periods}{tsp}(x)
    \method{periods}{tframe}(x)
    tframePeriods(x)
    \method{tframePeriods}{default}(x)
    \method{tframePeriods}{tftframe}(x)
    \method{tframePeriods}{tstframe}(x)
    \method{tframePeriods}{stamped}(x)
}

\arguments{
    \item{x}{
        A tframe or a tframed object.
    }

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
}
\examples{
    z <- ts(rnorm(100), start=c(1982,1), frequency=12)
    start(z)
    end(z)
    periods(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}


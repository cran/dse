\name{mineStepwise}
\alias{mineStepwise}
\title{Mine Stepwise}
\description{mineStepwise}
\usage{
    mineStepwise(data, essential.data=1,
      method="efroymson", f.crit=2, intercept=TRUE,
      subtract.means=FALSE,  standardize=FALSE, 
      lags.in=6, lags.out=6, trend=FALSE, plot.=TRUE)
}
\arguments{
    \item{data}{TSdata}
    \item{essential.data}{An integer vector indication important data.}
    \item{method}{method to pass to stepwise.}
    \item{f.crit}{See details.}
    \item{intercept}{See details.}
    \item{subtract.means}{See details.}
    \item{standardize}{See details.}
    \item{lags.in}{See details.}
    \item{lags.out}{See details.}
    \item{trend}{See details.}
    \item{plot.}{See details.}
}
\value{x}
\details{This documentation is out of date.
    Data should be of class TSdata. essential.data must have length 1.
    standardize and subtract means ...
    The result is a list with the results of stepwise,..., and
    several vectors indicating information about the columns of the matrix
    passed to stepwise: 
    io.indicator indicating an input (FALSE) or output (TRUE) variable
    v.indicator  indicating which series
    lag.indicator indicating the lag
    s.input.indicator and s.output.indicator are logic matrices
    \code{length(stepwise$rss)} by m and p respectively indicating if
    a series is included for each element of rss.
}

\concept{DSE}
\keyword{ts}


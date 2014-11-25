\name{tfwindow}
\alias{tfwindow}
\alias{tfwindow.default}
\alias{tfwindow.ts}

\title{Truncate a Time Series}
\description{Truncate a time series object to a time window.}
\usage{
    tfwindow(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
    \method{tfwindow}{default}(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
    \method{tfwindow}{ts}(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
    }

\arguments{
    \item{x}{A time series object.}
    \item{start.}{A start date of a format compatible with the time series}
    \item{end.}{An end date of a format compatible with the time series}
    \item{tf}{A tframe or tframed object}
    \item{warn}{A logical indicating if warning should be produced}
}
\value{
A time series object similar to x, but typically spanning a shorter time period.
}
\details{
If \code{start.} or \code{end.} are omitted and \code{tf} is specified then the
start or end is taken from the \code{tf} object.
For ts class objects this function calls window but makes more effort to 
preserve \code{seriesNames} if \code{x} has them. 
It also supports the optional argument \code{warn} to suppress
warning messages. Frequently it is convenient to write code which always 
truncates to a window without first checking if the data is already within 
the window. Since window produces a warning in this situation, the optional 
argument is frequently useful when tfwindow is used by other code.
In Splus \code{tfwindow} also corrects for some bugs in older versions 
of code{window}.
}
\examples{
  z <- ts(matrix(rnorm(24),24,1), start=c(1980,1), frequency=4)
  zz <- tfwindow(z, start=c(1982,2))
  zzz <- matrix(rnorm(24),24,1)
  tframe(zzz) <- tframe(z)
  tfwindow(zzz, tf=tframe(zz))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



\name{earliest.end}
\alias{earliest.end}
\alias{earliest.end.index}
\alias{earliest.end.index.tframe}
\alias{earliest.end.index.tframe.default}
\alias{earliest.end.index.tframe.tstframe}
\alias{earliest.end.tframe}
\alias{earliest.start}
\alias{earliest.start.index}
\alias{earliest.start.index.tframe}
\alias{earliest.start.index.tframe.default}
\alias{earliest.start.index.tframe.tstframe}
\alias{earliest.start.tframe}
\alias{latest.end}
\alias{latest.end.index}
\alias{latest.end.index.tframe}
\alias{latest.end.index.tframe.default}
\alias{latest.end.index.tframe.tstframe}
\alias{latest.end.tframe}
\alias{latest.start}
\alias{latest.start.index}
\alias{latest.start.index.tframe}
\alias{latest.start.index.tframe.default}
\alias{latest.start.index.tframe.tstframe}
\alias{latest.start.tframe}

\title{Start and End for Multiple Time Series Objects}
\description{Return start or end date (or index of the object) from multiple
    time series objects.
}

\usage{
    earliest.end(x, ...)
    earliest.end.tframe(x, ...)
    earliest.end.index(x, ...)
    earliest.end.index.tframe(x, ...)
    \method{earliest.end.index.tframe}{default}(x, ...)
    \method{earliest.end.index.tframe}{tstframe}(x, ...)

    earliest.start(x, ...)
    earliest.start.tframe(x, ...)
    earliest.start.index(x, ...)
    earliest.start.index.tframe(x, ...)
    \method{earliest.start.index.tframe}{default}(x, ...)
    \method{earliest.start.index.tframe}{tstframe}(x, ...)

    latest.end(x, ...)
    latest.end.tframe(x, ...)
    latest.end.index(x, ...)
    latest.end.index.tframe(x, ...)
    \method{latest.end.index.tframe}{default}(x, ...)
    \method{latest.end.index.tframe}{tstframe}(x, ...)

    latest.start(x, ...)
    latest.start.tframe(x, ...)
    latest.start.index(x, ...)
    latest.start.index.tframe(x, ...)
    \method{latest.start.index.tframe}{default}(x, ...)
    \method{latest.start.index.tframe}{tstframe}(x, ...)

}

\arguments{
    \item{x}{A tframe or tframed object.}
    \item{...}{Additional tframe or tframed objects.}
}
\value{A date or index.}
\details{
    These functions calculate the start
    and end of each object in the argument and return a result by 
    comparing across objects. Thus, latest.start returns the start date of
    the object which starts latest and latest.start.index returns the 
    corresponding index of the object in the argument list. 
}
\seealso{
    \code{\link{tframe}}
    \code{\link{tfwindow}}
    \code{\link{tftruncate}}
    \code{\link{trim.na}}
}

\examples{
    t1<-ts(c(1,2,3,4,5), start=c(1991,1))
    t2<-ts(c(2,3,4,5,6,7,8), start=c(1992,1))
    t3<-ts(c(NA,2,3,4,5), start=c(1991,1))

    latest.start(t1,t2,t3)  # 1992 1 corresponding to the starting date of 
                            # the object which starts latest (t2)
    latest.start(t1,t3)     # both start in 1991 1 (NAs count as data)
    latest.start(tbind(t1,t2,t3)) # tbind gives a single object starting in 1991 1
    latest.start(t2, tbind(t1,t2,t3))

    latest.start.index(t1,t2,t3)  # position of t2 in the argument list
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}


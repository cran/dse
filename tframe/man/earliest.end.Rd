\name{earliest.end}
\alias{earliest.end}
\alias{earliest.endIndex}
\alias{tframeEarliestEndIndex}
\alias{tframeEarliestEndIndex.default}
\alias{tframeEarliestEndIndex.tstframe}
\alias{earliest.tframeEnd}
\alias{earliest.start}
\alias{earliest.startIndex}
\alias{tframeEarliestStartIndex}
\alias{tframeEarliestStartIndex.default}
\alias{tframeEarliestStartIndex.tstframe}
\alias{earliest.tframeStart}
\alias{latest.end}
\alias{latest.endIndex}
\alias{tframeLatestEndIndex}
\alias{tframeLatestEndIndex.default}
\alias{tframeLatestEndIndex.tstframe}
\alias{latest.tframeEnd}
\alias{latest.start}
\alias{latest.startIndex}
\alias{tframeLatestStartIndex}
\alias{tframeLatestStartIndex.default}
\alias{tframeLatestStartIndex.tstframe}
\alias{latest.tframeStart}

\title{Start and End for Objects with Multiple Time Series}
\description{Return start or end date (or index of the object) from multiple
    time series objects.
}

\usage{
    earliest.end(x, ...)
    earliest.tframeEnd(x, ...)
    earliest.endIndex(x, ...)
    tframeEarliestEndIndex(x, ...)
    \method{tframeEarliestEndIndex}{default}(x, ...)
    \method{tframeEarliestEndIndex}{tstframe}(x, ...)

    earliest.start(x, ...)
    earliest.tframeStart(x, ...)
    earliest.startIndex(x, ...)
    tframeEarliestStartIndex(x, ...)
    \method{tframeEarliestStartIndex}{default}(x, ...)
    \method{tframeEarliestStartIndex}{tstframe}(x, ...)

    latest.end(x, ...)
    latest.tframeEnd(x, ...)
    latest.endIndex(x, ...)
    tframeLatestEndIndex(x, ...)
    \method{tframeLatestEndIndex}{default}(x, ...)
    \method{tframeLatestEndIndex}{tstframe}(x, ...)

    latest.start(x, ...)
    latest.tframeStart(x, ...)
    latest.startIndex(x, ...)
    tframeLatestStartIndex(x, ...)
    \method{tframeLatestStartIndex}{default}(x, ...)
    \method{tframeLatestStartIndex}{tstframe}(x, ...)

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
    the object which starts latest and latest.startIndex returns the 
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

    latest.startIndex(t1,t2,t3)  # position of t2 in the argument list
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}


\name{tfexpand}
\alias{tfexpand}
\alias{tfexpand.default}
\alias{tfexpand.tframe}
\alias{tframeExpand}
\alias{tframeExpand.default}
\alias{tframeExpand.tstframe}
\alias{tftruncate}
\alias{tftruncate.default}
\alias{tftruncate.tframe}
\alias{tframeTruncate}
\alias{tframeTruncate.default}
\alias{tframeTruncate.tstframe}

\title{Expand a Tframe or Tframed Object.}
\description{Expand a tframe or tframed object.}

\usage{
    tfexpand(x, add.start = 0, add.end = 0)
    \method{tfexpand}{default}(x, add.start = 0, add.end = 0)
    \method{tfexpand}{tframe}(x, add.start = 0, add.end = 0)
    tframeExpand(x, add.start=0, add.end=0)
    \method{tframeExpand}{default}(x, add.start=0, add.end=0)
    \method{tframeExpand}{tstframe}(x, add.start=0, add.end=0)

    tftruncate(x, start=NULL, end=NULL)
    \method{tftruncate}{default}(x, start=NULL, end=NULL)
    \method{tftruncate}{tframe}(x, start=NULL, end=NULL)
    tframeTruncate(x, start=NULL, end=NULL)
    \method{tframeTruncate}{default}(x, start=NULL, end=NULL)
    \method{tframeTruncate}{tstframe}(x, start=NULL, end=NULL)
}
\arguments{
    \item{x}{A tframe or tframed object.}
    \item{start}{an integer indicating the position at which the new tframe is
        to start.}
    \item{end}{an integer indicating the position at which the new tframe is
        to end.}
    \item{add.start}{an integer indicating the number of periods on the
        beginning.}
    \item{add.end}{an integer indicating the number of periods on the end.}
 }
\value{A tframe or tframed object.}
\details{
    These functions are like tfwindow but use position indicators (rather 
    than dates) and work with the tframe or tframed data. They are low level
    utilities for other functions.
}
\seealso{
    \code{\link{tfwindow}}
    \code{\link{tframed}}
}

\examples{
    tframe(tftruncate(ts(rnorm(100), start=c(1982,1), frequency=12), start=5))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



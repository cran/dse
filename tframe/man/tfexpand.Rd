\name{tfexpand}
\alias{tfexpand}
\alias{tfexpand.default}
\alias{tfexpand.tframe}
\alias{tfexpand.tframe.default}
\alias{tfexpand.tframe.tstframe}
\alias{tftruncate}
\alias{tftruncate.default}
\alias{tftruncate.tframe}
\alias{tftruncate.tframe.default}
\alias{tftruncate.tframe.tstframe}

\title{Expand a Tframe or Tframed Object.}
\description{Expand a tframe or tframed object.}

\usage{
    tfexpand(x, add.start=0, add.end=0)
    \method{tfexpand}{default}(x, start, end)
    tfexpand.tframe(x, add.start=0, add.end=0)
    \method{tfexpand.tframe}{default}(x, add.start=0, add.end=0)
    \method{tfexpand.tframe}{tstframe}(x, add.start=0, add.end=0)

    tftruncate(x, start=NULL, end=NULL)
    \method{tftruncate}{default}(x, start=NULL, end=NULL)
    tftruncate.tframe(x, start=NULL, end=NULL)
    \method{tftruncate.tframe}{default}(x, start=NULL, end=NULL)
    \method{tftruncate.tframe}{tstframe}(x, start=NULL, end=NULL)
}
\arguments{
    \item{x}{A tframe or tframed object.}
    \item{start}{
        An integer indicating the position at which the new tframe is
        to start.
    }
    \item{end}{
        An integer indicating the position at which the new tframe is
        to end.
    }
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



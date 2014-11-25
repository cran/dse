\name{periods.TSdata}
\alias{periods.TSdata}
\alias{periods.TSestModel}
\alias{start.TSdata}
\alias{end.TSdata}
\alias{frequency.TSdata}
\alias{tframed.TSdata}
\alias{tfwindow.TSdata}
\alias{tbind.TSdata}
\alias{trim.na.TSdata}

\alias{window.TSdata}

\title{Specific Methods for tframed Data}
\description{See the generic function description.}
\usage{
    \method{periods}{TSdata}(x)
    \method{periods}{TSestModel}(x)
    \method{start}{TSdata}(x)
    \method{end}{TSdata}(x)
    \method{frequency}{TSdata}(x)
    \method{tframed}{TSdata}(x, tf=NULL, input.names=NULL, output.names=NULL)
    \method{tfwindow}{TSdata}(x, start=NULL, end=NULL, tf=NULL, warn=T)
    \method{tbind}{TSdata}(d1, d2)
    \method{trim.na}{TSdata}(x, start.=T, end.=T)
}

\seealso{
    \code{\link{periods}}, 
    \code{\link{start}},
    \code{\link{end}},
    \code{\link{frequency}},
    \code{\link{tframed}},
    \code{\link{tfwindow}},
    \code{\link{tbind}},
    \code{\link{trim.na}}
}

%\keyword{DSE}
\keyword{ts}


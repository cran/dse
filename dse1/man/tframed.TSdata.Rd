\name{tframed.TSdata}
\alias{tframed.TSdata}
\alias{tframe<-.TSdata}
\alias{tfwindow.TSdata}
\alias{tbind.TSdata}
\alias{trim.na.TSdata}

\alias{window.TSdata}

\title{Specific Methods for tframed Data}
\description{See the generic function description.}
\usage{
    \method{tframed}{TSdata}(x, tf=NULL, input.names=NULL, output.names=NULL)
%    \method{tframe<-}{TSdata}(x, value)
    \method{tfwindow}{TSdata}(x, start=NULL, end=NULL, tf=NULL, warn=TRUE)
    \method{tbind}{TSdata}(d1, d2)
    \method{trim.na}{TSdata}(x, start.=TRUE, end.=TRUE)
}
\arguments{
    \item{x}{See the generic function.}
    \item{tf}{a time frame. See the generic function.}
    \item{input.names}{See the generic function.}
    \item{output.names}{See the generic function.}
    \item{start}{See the generic function.}
    \item{start.}{See the generic function.}
    \item{end}{See the generic function.}
    \item{end.}{See the generic function.}
    \item{d1}{See the generic function.}
    \item{d2}{See the generic function.}
    \item{warn}{logical indicating if some warning messages should be suppressed.}
}

\seealso{
    \code{\link{tframed}},
    \code{\link{tfwindow}},
    \code{\link{tbind}},
    \code{\link{trim.na}}
}

%\keyword{DSE}
\keyword{ts}


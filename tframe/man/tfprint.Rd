\name{tfprint}
\alias{tfprint}
\alias{tfprint.default}
\alias{tfprint.tframe}
\alias{tfprint.tframe.default}

\title{Print Tframed Objects}
\description{
    Print tframe or tframed objects. 
}
\usage{
    tfprint(x, ...)
    \method{tfprint}{default}(x, ...)
    \method{tfprint}{tframe}(x, ...)
    \method{tfprint}{tframe.default}(x, digits=NULL, quote=T, prefix="", ...)
}

\arguments{
    \item{x}{A tframe or tframed object.}
}
\value{tfprint methods return the object invisibly.}
\details{
    \code{tfprint} prints data in a tframed object while 
    \code{tfprint.tframe} prints the tframe.
    In many cases these are the same as print methods. 
    However, tfprint tries to provide an alternate generic
    mechanism that is consistent with the tframe view of the data.
    This may not always be the preferred print method. 
    Also, new classes of time series may define
    there own print methods in ways which use a different logic
    from the tframe library. Thus tfprint provides a
    way to program functions which use methods consistent with the tframe
    library logic.
}
\sideffects{An object is printed.}
\seealso{
    \code{\link{tfplot}}
    \code{\link{tframe}}
    \code{\link{tframed}}
    \code{\link{print}}
    \code{\link{plot}}
}
\examples{
    tfprint(ts(rnorm(100)))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}



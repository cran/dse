\name{tfprint}
\alias{tfprint}
\alias{tfprint.default}
\alias{tfprint.tframe}
\alias{tframePrint}
\alias{tframePrint.default}

\title{Print Tframed Objects}
\description{
    Print tframe or tframed objects. 
}
\usage{
    tfprint(x, ...)
    tframePrint(x, ...)
    \method{tfprint}{default}(x, ...)
    \method{tfprint}{tframe}(x, ...)
    \method{tframePrint}{default}(x, digits=NULL, quote=TRUE, prefix="", ...)
}

\arguments{
    \item{x}{a tframe or tframed object.}
    \item{digits}{see the generic method.}
    \item{quote}{see the generic method.}
    \item{prefix}{see the generic method.}
    \item{...}{arguments to be passed to other methods.}
}
\value{tfprint methods return the object invisibly.}
\details{
    \code{tfprint} prints data in a tframed object while 
    \code{tframePrint} prints the tframe.
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



\name{tframe}
\alias{tframe}
\alias{tframe.default}
\alias{tframe.ts}
\alias{tframe<-}
\alias{tframe<-.default}
\alias{tframe<-.rts}
\alias{tframe<-.cts}
\alias{tframe<-.its}
\alias{tframe<-.ts}
\alias{is.tframe}

\title{Extract or Set a tframe}
\description{
    The first usage extracts the tframe from a tframed object. 
    These are generic functions. The default method attempts to build a
    tframe for time series, matrix, or array which is not tframed.
    The second usage assigns the tframe to an object.

 There should be 
1/ methods for classes like TSdata, ts, rts
2/ the default should work for "true tframe" methods
3/ the default should kludge for old tsp


}

\usage{
    tframe(x)
    \method{tframe}{default}(x)
    \method{tframe}{ts}(x)

    tframe(x) <- value
    \method{tframe}{default}(x) <- value
    \method{tframe}{rts}(x) <- value
    \method{tframe}{cts}(x) <- value
    \method{tframe}{its}(x) <- value
    \method{tframe}{ts}(x) <- value

    is.tframe(tf)
}

\arguments{
    \item{x}{a tframed object or an object to be tframed.}
    \item{value}{a tframe.}
    \item{tf}{any object.}
}
\value{Depends.}
\details{
    The first usage returns the tframe of a tframed object. The second
    usage sets the tframe of an object. \code{is.tframe} returns a logical.
    
    The true tframe approach is to set a "tframe" attribute for an object. This
    attribute has a class which indicates the time framing which is used. The 
    the time frame information is often secondary, in the sense that it does not
    describe the object structure, but only provide some additional information
    which is useful for plots and when printing the object. By putting this in 
    an attribute, the objects class can be used in the more usual way, for indicating
    information about the structure of the object.  For
    these true tframe objects the default \code{tframe} and code{tframe<-} will
    often be adequate. The generic/method approach allows for special case (like
    TSdata where the tframe information is not an attribute of the object, but
    rather an attribute of the data matrices which are elements of the object).

    The generic/method approach also allows for (faking) tframe assignment and 
    extraction with classes like rts, ctc, its, ts, and others which may appear,
    that try to make the time description part of the object class. (Not 
    a "tframe" approach.) The problem is to extract real tframes and
    also fake these other classes and old style tsp objects so they look 
    like tfamed objects. Another approach would be to mutilate these objects 
    and force them really be tframed objects (to have a tframe attribute), 
    but that risks conflicting with other (non tframe) code which used 
    the objects. This faking is accomplished by specific methods of the classes,
    and for old style tsp objects it is built into the default.
}
\seealso{
    \code{\link{tframed}}
}
\examples{
    z <- tframe(ts(rnorm(100), start=c(1982,1), frequency=12))
    is.tframe(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



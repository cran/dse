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
\alias{tframed}
\alias{tframed.default}
\alias{is.tframe}
\alias{is.tframed}

\title{Extract or Set a tframe}
\description{Extract or set the tframe of an object.}

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

    tframed(x, tf=NULL, names = NULL)
    \method{tframed}{default}(x, tf = NULL, names = NULL)

    is.tframe(x)
    is.tframed(x)
}

\arguments{
    \item{x}{an object (to which a tframe is assigned in assignment methods).}
    \item{value}{a tframe.}
    \item{tf}{a tframe object or a tframed object from which a tframe is taken.}
    \item{names}{optional vector of strings to specify new series names.}
}
\value{Depends.}
\details{
    The first usage returns the tframe of a tframed object. The assignment
    methods and tframed set the tframe of an object. \code{is.tframe} returns 
    a logical.
    
    The pure tframe approach is to set a "tframe" attribute for an object. This
    attribute has a class which indicates the time framing which is used. The 
    the time frame information is often secondary, in the sense that it does not
    describe the object structure, but only provides some additional information
    which is useful for doing time based operations on the data, plotting, and
    printing the object. By putting this in 
    an attribute, the objects class can be used for indicating other
    information about the structure of the object.  For
    these pure tframe objects the default \code{tframe} and code{tframe<-} will
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

    This \code{tframed} constructor is simply a shortcut for assigning the tframe
    (tframe(x) <- tf) and series names (seriesNames(x) <-  names) to an object,
    but never assigns NULL values, so the result is guaranteed to be a 
    code{tframed} object. It is like \code{ts} but enables the tframe library's  
    methods for handling time. If the \code{tf} argument is a 
    \code{tframed} object rather than a \code{tframe}, then the code{tframe} 
    is extracted and used. If the \code{names} argument is not mode "character"
    of appropriate length, then \code{seriesNames(names)} is used. These make it
    simple to assign the time frame and names of one object to another by
    \code{z <- tframed(x, tf=y, names=y)}.
    
    \code{is.tframed} returns TRUE if a
    \code{tframe()} can extract a tframe from the object. This is true for many
    objects which are not truly tframed (like ts objects), since \code{tframe()}
    tries fairly hard to build a tframe for the object. 
}
\seealso{
    \code{\link{tfstart}},
    \code{\link{tfend}},
    \code{\link{tffrequency}},
    \code{\link{tfperiods}},
    \code{\link{tftime}},
    \code{\link{tfdiff}}
}
\examples{
    z <- tframe(ts(rnorm(100), start=c(1982,1), frequency=12))
    is.tframe(z)
    zz <- tframed(matrix(rnorm(200), 100,2), tf=z)
    is.tframed(zz)
    zzz <- tframed(matrix(rnorm(200), 100,2), tf=zz)
    is.tframed(zzz)
    tframe(zzz)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}


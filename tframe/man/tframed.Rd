\name{tframed}
\alias{tframed}
\alias{tframed.default}
\alias{is.tframed}
\alias{is.Ttframed}

\title{Construct a Tframed Object}
\description{
    Both functions construct a tframed object. 
    These are generic functions. If tf is not a tframe object but has the 
    necessary information then the default methods attempt to build 
    a tframe.
}

\usage{
    tframed(x, tf=NULL, names = NULL)
    \method{tframed}{default}(x, tf = NULL, names = NULL)

    is.tframed(x)
    is.Ttframed(x)
}

\arguments{
    \item{x}{a tframed object or an object to be tframed.}
    \item{tf}{a tframe.}
    \item{names}{optional vector of strings to specify new series names.}
    \item{value}{a tframe attribute to be applied to x.}
}
\value{A tframed object.}
\details{
    This \code{tframed}constructor is like ts but enables the tframe library's  
    generic methods for handling time. \code{is.tframed} returns TRUE if a
    \code{tframe()} can extract a tframe from the object. This is true for many
    objects which are not truly tframed (like ts objects), since \code{tframe()}
    tries fairly hard to build a tframe for the object. \code{is.Ttframed} 
    returns TRUE only if the object is truly tframed (i.e. has an attribute
    \code{tframe}).
}
\seealso{
    \code{\link{tframe}}
}
\examples{
    z <- tframed(matrix(rnorm(200), 100,2), 
            tf=list(start=c(1982,1), frequency=12))
    is.tframed(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



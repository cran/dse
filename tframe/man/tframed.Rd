\name{tframed}
\alias{tframed}
\alias{tframed.default}
\alias{settf}
\alias{settf.default}
\alias{settf.tftframe}
\alias{is.tframed}

\title{Construct a Tframed Object}
\description{
    Both functions construct a tframed object. 
    These are generic functions. If tf is not a tframe object but has the 
    necessary information then the default methods attempt to build 
    a tframe.
}

\usage{
    tframed(x, ...)
    \method{tframed}{default}(x, tf = NULL, names = NULL)

    settf(value, x)
    \method{settf}{default}(value, x)
    \method{settf}{tftframe}(value, x)

    is.tframed(x)
}

\arguments{
    \item{x}{a tframed object or an object to be tframed.}
    \item{tf}{a tframe.}
    \item{names}{optional vector of strings to specify new series names.}
    \item{value}{a tframe attribute to be applied to x.}
    \item{...}{arguments to be passed to other methods.}
}
\value{A tframed object.}
\details{
    This constructor is like ts but enables the tframe library's generic 
    methods for handling time.
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



\name{lag}
\alias{lag}
\alias{lag.tframe}
\alias{tsmatrix}
\alias{tsmatrix.default}
\title{Misc Generic Replacements}
\description{
    These functions replace functions in S or R with generics if necessary.}
\usage{
    lag(x, ...)
    \method{lag}{default}(x, ...)
    \method{lag}{tframe}(x, ...)

    tsmatrix(x, ...)
    \method{tsmatrix}{default}(x, ...)
}
\details{
    The default for the generic is set to the original function. See the S/R
    documentation for the original.
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



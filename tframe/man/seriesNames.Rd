\name{seriesNames}
\alias{seriesNames}
\alias{seriesNames.default}
\alias{seriesNames<-}
\alias{seriesNames<-.default}

\title{Names of Series in a time series object}
\description{Extract or set names of series in a time series object.}
\usage{
    seriesNames(x)
    \method{seriesNames}{default}(x)

    seriesNames(data) <- names 
}
\arguments{
    \item{x}{a time series object.}
    \item{data}{a time series object.}
}
\value{
The first usage returns a vector of strings with the series names. 
The assignment method makes \code{names} (a vector of strings) 
the series names of data.
}
\seealso{
    \code{\link{tframed}},
    \code{\link[dse1]{seriesNamesInput}},
    \code{\link[dse1]{seriesNamesOutput}}
}

\examples{
    z <- matrix(rnorm(100), 50,2)
    seriesNames(z) <- c("a", "b")
    seriesNames(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}



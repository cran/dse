\name{series.names}
\alias{series.names}
\alias{series.names.default}
\alias{series.names<-}
\alias{series.names<-.default}

\title{Names of Series in a time series object}
\description{Extract or set names of series in a time series object.}
\usage{
    series.names(x)
    \method{series.names}{default}(x)

    series.names(data) <- names 
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
    \code{\link[dse1]{input.series.names}},
    \code{\link[dse1]{output.series.names}}
}

\examples{
    z <- matrix(rnorm(100), 50,2)
    series.names(z) <- c("a", "b")
    series.names(z)
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}



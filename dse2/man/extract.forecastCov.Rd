\name{extract.forecastCov}
\alias{extract.forecastCov}
\alias{extract.forecastCov.estimatorsWRTdata}
\alias{extract.forecastCov.estimators.from.model}

\title{Extract Forecast Covariance}
\description{extract.forecastCov}
\usage{
    extract.forecastCov(e, n)
    \method{extract.forecastCov}{estimatorsWRTdata}(e, n)
    \method{extract.forecastCov}{estimators.from.model}(e, n)
}
\arguments{
    \item{e}{A forecastCov object.}
    \item{n}{A vector on integers.}
}
\value{A forecastCov object.}
\details{
    Select a subset of models and their forecast covariances from a
    larger object.
}
\seealso{
    \code{\link{forecastCov}}
}

%\keyword{DSE}
\keyword{ts}


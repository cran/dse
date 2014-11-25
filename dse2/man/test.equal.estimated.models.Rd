\name{test.equal.estimated.models}
\alias{test.equal.estimated.models}
\alias{test.equal.EstEval}
\alias{test.equal.forecast}
\alias{test.equal.forecastCov}
\alias{test.equal.horizonForecasts}
\alias{test.equal.MonteCarloSimulations}

\title{Specific Methods for Testing Equality}
\description{See the generic function description.}
\usage{
    \method{test.equal}{estimated.models}(obj1, obj2, fuzz = 0)
    \method{test.equal}{EstEval}(obj1, obj2, fuzz=0)
    \method{test.equal}{forecast}(obj1, obj2, fuzz=1e-14)
    \method{test.equal}{forecastCov}(obj1, obj2, fuzz=1e-14)
    \method{test.equal}{horizonForecasts}(obj1, obj2, fuzz=1e-14)
    \method{test.equal}{MonteCarloSimulations}(obj1, obj2, fuzz=1e-16)
}
\arguments{
    \item{obj1}{an object which is to be compared with the second object.}
    \item{obj2}{an object which is to be compared with the first object.}
    \item{fuzz}{tolerance for numerical comparisons. Values within fuzz will be
       considered equal.}
}

\seealso{ \code{\link{test.equal}} }

%\keyword{DSE}
\keyword{ts}


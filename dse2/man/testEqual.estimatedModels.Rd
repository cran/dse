\name{testEqual.estimatedModels}
\alias{testEqual.estimatedModels}
\alias{testEqual.EstEval}
\alias{testEqual.forecast}
\alias{testEqual.forecastCov}
\alias{testEqual.horizonForecasts}
\alias{testEqual.MonteCarloSimulations}

\title{Specific Methods for Testing Equality}
\description{See the generic function description.}
\usage{
    \method{testEqual}{estimatedModels}(obj1, obj2, fuzz = 0)
    \method{testEqual}{EstEval}(obj1, obj2, fuzz=0)
    \method{testEqual}{forecast}(obj1, obj2, fuzz=1e-14)
    \method{testEqual}{forecastCov}(obj1, obj2, fuzz=1e-14)
    \method{testEqual}{horizonForecasts}(obj1, obj2, fuzz=1e-14)
    \method{testEqual}{MonteCarloSimulations}(obj1, obj2, fuzz=1e-16)
}
\arguments{
    \item{obj1}{an object which is to be compared with the second object.}
    \item{obj2}{an object which is to be compared with the first object.}
    \item{fuzz}{tolerance for numerical comparisons. Values within fuzz will be
       considered equal.}
}

\seealso{ \code{\link[tframe]{testEqual}} }

\concept{DSE}
\keyword{ts}


\name{TSdata.forecastCov}
\alias{TSdata.forecastCov}
\alias{TSestModel.coefEstEval}
\alias{TSmodel.forecastCov}
\alias{TSmodel.coefEstEval}

\title{TS Extractor Specific Methods}
\description{See the generic function description.}
\usage{
    \method{TSdata}{forecastCov}(data, ...)
    \method{TSestModel}{coefEstEval}(obj)
    \method{TSmodel}{forecastCov}(obj, select=1, ...)
    \method{TSmodel}{coefEstEval}(obj, ...)
}
\arguments{
    \item{data}{an object from which to extract the TSdata.}
    \item{obj}{an object from which to extract the TSmodel or TSestModel.}
    \item{select}{an integer indicating which of multiple models to extract.}
    \item{...}{arguments to be passed to other methods.}
}

\seealso{ 
    \code{\link[dse1]{TSdata}} 
    \code{\link[dse1]{TSestModel}} 
    \code{\link[dse1]{TSmodel}} 
}

%\keyword{DSE}
\keyword{ts}


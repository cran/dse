\name{TSdata.forecastCov}
\alias{TSdata.forecastCov}
\alias{TSestModel.coef.ee}
\alias{TSmodel.forecastCov}
\alias{TSmodel.coef.ee}

\title{TS Extractor Specific Methods}
\description{See the generic function description.}
\usage{
    \method{TSdata}{forecastCov}(obj)
    \method{TSestModel}{coef.ee}(obj)
    \method{TSmodel}{forecastCov}(obj, select=1)
    \method{TSmodel}{coef.ee}(obj)
}
\arguments{
    \item{obj}{an object from which to extract the TSdata, TSmodel, or
       TSestModel.}
    \item{select}{an integer indicating which of multiple models to extract.}
}

\seealso{ 
    \code{\link{TSdata}} 
    \code{\link{TSestModel}} 
    \code{\link{TSmodel}} 
}

%\keyword{DSE}
\keyword{ts}


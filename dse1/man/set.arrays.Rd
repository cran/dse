\name{set.arrays}
\alias{set.arrays}
\alias{set.arrays.ARMA}
\alias{set.arrays.SS}
\alias{set.arrays.TSestModel}

\title{Set TSmodel Array Information}
\description{
    Complete parameter array information based on parameter vector settings.
    This function is used internally and is not normally called by a user. 
}
\usage{
    set.arrays(model, parms=NULL)
    \method{set.arrays}{ARMA}(model, parms=NULL)
    \method{set.arrays}{SS}(model, parms=NULL)
    \method{set.arrays}{TSestModel}(model, parms=NULL)
    }
\arguments{
    \item{model}{An object of class TSmodel.}
}
\value{A TSmodel object.}

%\keyword{DSE}
\keyword{ts}


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
    set.arrays(model, coefficients=NULL)
    \method{set.arrays}{ARMA}(model, coefficients=NULL)
    \method{set.arrays}{SS}(model, coefficients=NULL)
    \method{set.arrays}{TSestModel}(model, coefficients=NULL)
    }
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{coefficients}{A vector of new values for the model coefficients
    (parameters).}
}
\value{A TSmodel object.}

%\keyword{DSE}
\keyword{internal}



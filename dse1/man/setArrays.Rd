\name{setArrays}
\alias{setArrays}
\alias{setArrays.ARMA}
\alias{setArrays.SS}
\alias{setArrays.TSestModel}
\alias{setArrays.TSrestrictedModel}

\title{Set TSmodel Array Information}
\description{
    Complete parameter array information based on parameter vector settings.
    This function is used internally and is not normally called by a user. 
}
\usage{
    setArrays(model, coefficients=NULL)
    \method{setArrays}{ARMA}(model, coefficients=NULL)
    \method{setArrays}{SS}(model, coefficients=NULL)
    \method{setArrays}{TSestModel}(model, coefficients=NULL)
    \method{setArrays}{TSrestrictedModel}(model, coefficients=coef(model))
    }
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{coefficients}{A vector of new values for the model coefficients
    (parameters).}
}
\value{A TSmodel object.}

\concept{DSE}
\keyword{internal}



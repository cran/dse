\name{coef.TSmodel}
\alias{coef.TSmodel}
\alias{coef.TSestModel}
\alias{coef.TSrestrictedModel}

\alias{coef<-}
\alias{coef<-.default}

\title{Extract Model Parameters}
\description{
    Set or extract coefficients (parameter values) of model objects.
}
\usage{
    \method{coef}{TSmodel}(object, ...)
    \method{coef}{TSestModel}(object, ...)
    \method{coef}{TSrestrictedModel}(object, ...)
    coef(object) <- value
    \method{coef}{default}(object) <- value
    }
\arguments{
    \item{object}{An object of class TSmodel or TSestModel.}
    \item{value}{value to be assigned to object.}
    \item{...}{(further arguments, currently disregarded).}
    }
\value{A vector of parameter values.}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- estVARXls(eg1.DSE.data.diff)
    coef(model)
    coef(model) <- 0.1 + coef(model)
}
\concept{DSE}
\keyword{ts}



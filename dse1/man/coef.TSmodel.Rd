\name{coef.TSmodel}
\alias{coef.TSmodel}
\alias{coef.TSestModel}

\alias{coef<-}
\alias{coef<-.default}

\title{Extract Model Parameters}
\description{
    Set or extract coefficients (parameter values) of model objects.
}
\usage{
    \method{coef}{TSmodel}(object, ...)
    \method{coef}{TSestModel}(object, ...)
    coef(object) <- value
    \method{coef<-}{default}(object) <- value
    }
\arguments{
    \item{object}{An object of class TSmodel or TSestModel.}
    \item{...}{(further arguments, currently disregarded).}
    }
\value{A vector of parameter values.}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    coef(model)
    coef(model) <- 0.1 + coef(model)
}
%\keyword{DSE}
\keyword{ts}



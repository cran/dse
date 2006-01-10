\name{setTSmodelParameters}
\alias{setTSmodelParameters}
\alias{setTSmodelParameters.default}
\alias{setTSmodelParameters.ARMA}
\alias{setTSmodelParameters.SS}
\alias{setTSmodelParameters.TSrestrictedModel}

\title{Set TSmodel Parameter Information}
\description{
    Complete parameter vector information based on parameter array settings.
    This function is used internally and is not normally called by a user. 
}
\usage{
    setTSmodelParameters(model, constants=model$constants)
    \method{setTSmodelParameters}{default}(model, constants=TSmodel(model)$constants)
    \method{setTSmodelParameters}{ARMA}(model, constants=model$constants)
    \method{setTSmodelParameters}{SS}(model, constants=model$constants)
    \method{setTSmodelParameters}{TSrestrictedModel}(model,
         constants=TSmodel(model$TSmodel)$constants)  
}
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{constants}{A list of logical arrays indicating TRUE for any model
    array entries that should be treated as constants.}
}
\value{An object of class 'TSmodel'.}
\seealso{

\code{\link{setArrays}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- estVARXls(eg1.DSE.data.diff)
    model <- setTSmodelParameters(model)
}
\concept{DSE}
\keyword{ts}


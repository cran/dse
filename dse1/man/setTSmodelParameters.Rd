\name{setTSmodelParameters}
\alias{setTSmodelParameters}
\alias{setTSmodelParameters.default}
\alias{setTSmodelParameters.ARMA}
\alias{setTSmodelParameters.SS}

\title{Set TSmodel Parameter Information}
\description{
    Complete parameter vector information based on parameter array settings.
    This function is used internally and is not normally called by a user. 
}
\usage{
    setTSmodelParameters(model, constants=NULL)
    \method{setTSmodelParameters}{default}(model, constants=NULL)
    \method{setTSmodelParameters}{ARMA}(model, constants=NULL)
    \method{setTSmodelParameters}{SS}(model, constants=NULL)
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
%\keyword{DSE}
\keyword{ts}


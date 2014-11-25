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
    setTSmodelParameters(model)
    \method{setTSmodelParameters}{default}(model)
    \method{setTSmodelParameters}{ARMA}(model)
    \method{setTSmodelParameters}{SS}(model)
}
\arguments{
    \item{model}{An object of class TSmodel.}
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


\name{to.SSChol}
\alias{to.SSChol}
\alias{to.SSChol.TSmodel}
\alias{to.SSChol.TSestModel}

\title{Convert to Non-Innovation State Space Model}
\description{

Convert to a non-innovations state space representation using 
the given matrix as the measurement noise covariance. 
This function may not be working properly.
}
\usage{
    to.SSChol(model, ...)
    \method{to.SSChol}{TSmodel}(model, Om=diag(1,output.dimension(model)))
    \method{to.SSChol}{TSestModel}(model, Om=NULL)
}
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{Om}{
      a matrix to be used as the measurement noise covariance. If Om is
      not supplied and model is of class TSestModel then 
      \code{model$estimates$cov} is used. Otherwise, Om is set to the 
      identity matrix.}
    \item{...}{arguments to be passed to other methods.}
}
\value{
An object of class 'SS' 'TSmodel' containing a state space model which is
not in innovations form.
}
\details{
Convert to a non-innovations SS  representation using a Cholesky 
decomposition of Om as the coefficient matrix of the output noise.
}
\seealso{
\code{\link{to.SS.innov}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    model <- to.SSChol(model)
}
%\keyword{DSE}
\keyword{ts}


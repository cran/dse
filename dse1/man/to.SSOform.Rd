\name{to.SSOform}
\alias{to.SSOform}
\alias{to.SSOform.TSmodel}
\alias{to.SSOform.TSestModel}

\title{Convert to Oform}
\description{Convert a state space model to (observability?) form.}
\usage{
    to.SSOform(model)
    \method{to.SSOform}{TSmodel}(model)
    \method{to.SSOform}{TSestModel}(model)
    }
\arguments{
    \item{model}{An object of class TSmodel.}
}
\value{
An object of class 'SS' 'TSmodel' containing a state space 
model in observability form (more or less).
 WARNING: This function does not work properly.}
\seealso{
    \code{\link{to.SS.innov}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    model <- to.SSOform(model)
}
%\keyword{DSE}
\keyword{ts}


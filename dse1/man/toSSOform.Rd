\name{toSSOform}
\alias{toSSOform}
\alias{toSSOform.TSmodel}
\alias{toSSOform.TSestModel}

\title{Convert to Oform}
\description{Convert a state space model to (observability?) form.}
\usage{
    toSSOform(model)
    \method{toSSOform}{TSmodel}(model)
    \method{toSSOform}{TSestModel}(model)
    }
\arguments{
    \item{model}{An object of class TSmodel.}
}
\value{
An object of class 'SS' 'TSmodel' containing a state space 
model in observability form (more or less).
 WARNING: This function does not work properly.}
\seealso{
    \code{\link{toSSinnov}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- estVARXls(eg1.DSE.data.diff)
%    model <- toSSOform(model)
}
%\keyword{DSE}
\keyword{ts}


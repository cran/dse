\name{stability}
\alias{stability}
\alias{stability.ARMA}
\alias{stability.roots}
\alias{stability.TSmodel}
\alias{stability.TSestModel}

\title{Calculate Stability of a TSmodel}
\description{Calculate roots and their modulus and indicate stability.}
\usage{
    stability(obj, fuzz=1e-4, digits=8, verbose=TRUE)
    \method{stability}{ARMA}(obj, fuzz=1e-4, digits=8, verbose=TRUE)
    \method{stability}{roots}(obj, fuzz=1e-4, digits=8, verbose=TRUE)
    \method{stability}{TSmodel}(obj, fuzz=1e-4, digits=8, verbose=TRUE)
    \method{stability}{TSestModel}(obj, fuzz=1e-4, digits=8, verbose=TRUE)
    }
\arguments{
    \item{obj}{An object of class TSmodel.}
    \item{fuzz}{Roots within fuzz are considered equal.}
    \item{digits}{Printing precision.}
    \item{verbose}{Print roots and there moduli.}
}
\value{TRUE or FALSE if the model is stable or not stable.}
\sideffects{
The eigenvalues of the state transition matrix or the roots of the
determinant of the AR polynomial are printed if verbose is T.
}
\seealso{
\code{\link{McMillanDegree}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- estVARXls(eg1.DSE.data.diff)
stability(model)
}
\concept{DSE}
\keyword{ts}


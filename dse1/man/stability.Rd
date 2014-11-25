\name{stability}
\alias{stability}
\alias{stability.ARMA}
\alias{stability.SS}
\alias{stability.roots}
\alias{stability.TSmodel}
\alias{stability.TSestModel}

\title{Calculate stability of a TSmodel}
\description{Calculate roots and their modulus and indicate stability.}
\usage{
	stability(model, fuzz=1e-4, digits=8, verbose=T)}
\arguments{
\item{ model}{An object of class TSmodel.}
\item{ fuzz}
{Roots within fuzz are considered equal.}
\item{ digits}{Printing precision.}
\item{ verbose}{Print roots and there moduli.}
}
\value{True or False if the model is stable or not stable.}
\sideffects{
The eigenvalues of the state transition matrix or the roots of the
determinant of the AR polynomial are printed if verbose is T.
}
\seealso{

\code{\link{McMillan.degree}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
stability(model)
}
\keyword{DSE}
\keyword{ts}


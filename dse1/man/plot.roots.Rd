\name{plot.roots}
\alias{plot.roots}
\title{Plot Model Roots}
\description{Calculate and plot roots of a model.}
\usage{
	plot.roots(v, fuzz=0)}
\arguments{
	\item{v}{An object of class TSmodel or a vector of complex (or real) 
values as returned by the function roots.}
	\item{fuzz}{If non-zero then roots within fuzz distance are considered equal.}
}
\sideffects{The roots and a unit circle are plotted on the complex plane.}
\value{
The eigenvalues of the state transition matrix or the inverse of the roots of the
determinant of the AR polynomial are returned invisibly.
}
\seealso{
	\code{\link{add.plot.roots}}
	\code{\link{roots}}
	\code{\link{stability}}
	\code{\link{McMillan.degree}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.VARX.ls(eg1.DSE.data.diff)
	plot.roots(model)
}
\keyword{DSE}
\keyword{ts}


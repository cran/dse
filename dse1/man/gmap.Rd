\name{gmap}
\alias{gmap}
\title{Basis transformation of a model.}
\description{
Transform the basis for the state by a given invertible matrix.
}
\usage{
	gmap(g, model)
}
\arguments{
\item{ g}
{An invertible matrix}
\item{ model}
{
An object of class TSmodel.}

}
\value{

An equivalent model transformed using g.
}
\details{

If the input model is in state space form g is a change of basis 
for the state. If the input model is in ARMA form then the 
polynomials are premultiplied by g. If g is a scalar it is treated 
as a diagonal matrix.
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
gmap(2, model)
}
\keyword{DSE}
\keyword{ts}


\name{to.parsim}
\alias{to.parsim}
\alias{to.parsim.TSmodel}
\alias{to.parsim.TSestModel}

\title{convert to parsimonious form}
\description{Try to find a more parsimonious parameterization.}
\usage{
	to.parsim(model, minf = nlmin, small=1e-5, equiv=F,max.iter=100,max.fcal=200)
}
\arguments{
\item{ model}{An object of class TSmodel.}
\item{ arg}
{This function attempts to convert a model to a more parsimonious representation.
It does not work very well.}
}
\value{A TSmodel.}
{WARNING: This function does not work.}
\seealso{
\code{\link{to.SS.Oform}}
\code{\link{to.SS.innov}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.ARMA(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
#model <- to.parsim(model) not working
}
\keyword{DSE}
\keyword{ts}


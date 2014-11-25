\name{to.ARMA}
\alias{to.ARMA}
\alias{to.ARMA.ARMA}
\alias{to.ARMA.SS}
\alias{to.ARMA.TSestModel}

\title{Convert to an ARMA Model}
\description{
Convert a state space model to an ARMA representation. The state is 
eliminated by a method which uses an equivalence that can be demonstrated
by the Cayley Hamilton theorem.
}
\usage{
	to.ARMA(model)
}
\arguments{
\item{ model}
{
An object of class TSmodel.}

}
\value{

An object of class 'ARMA' 'TSmodel' containing an ARMA model.
}
\section{ REFERENCES} 

See for example Aoki(1989)
}
\seealso{

\code{\link{to.SS}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
model <- to.ARMA(model)
}
\keyword{DSE}
\keyword{ts}


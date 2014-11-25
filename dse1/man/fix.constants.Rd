\name{fix.constants}
\alias{fix.constant}
\title{Fix TSmodel Parameters to Constants}
\description{
	Fix any parameters within fuzz of 0.0 or 1.0 to exactly 0.0 or 1.0.
	This will not change the model much but will affect some estimation 
	techniques and information criteria results, as these are considered to 
	be constants rather than parameters.
}
\usage{
	fix.constants(model, fuzz=1e-5)
}
\arguments{
	\item{model}{An object of class TSmodel.}
	\item{fuzz}{Absolute difference to be considered equivalent.}}
\value{
	An object of class 'SS' 'TSmodel' with parameters near 0.0 or 1.0 set 
	as constants 0.0 or 1.0.
}

\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- to.ARMA(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
	model <- fix.constants(model)
}
\seealso{
	\code{\link{fixF}}
}
\keyword{DSE}
\keyword{ts}



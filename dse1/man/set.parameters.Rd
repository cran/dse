\name{set.parameters}
\alias{set.parameters}
\alias{set.parameters.TSmodel}
\alias{set.parameters.TSmodel.ARMA}
\alias{set.parameters.TSmodel.SS}

\title{Set TSmodel Parameter Information}
\description{
	Complete parameter vector information based on parameter array settings.
	This function is used internally and is not normally called by a user. 
}
\usage{
	set.parameters(model)
}
\arguments{
	\item{ model}{An object of class TSmodel.}
}
\value{An object of class 'TSmodel'.}
\seealso{

\code{\link{set.arrays}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.VARX.ls(eg1.DSE.data.diff)
	model <- set.parameters(model)
}
\keyword{DSE}
\keyword{ts}


\name{reconstitute}
\alias{reconstitute}
\alias{reconstitute.default}
\alias{reconstitute.concentrate}
\alias{reconstitute.TSdataconcentrate}
\alias{is.TSdatareconstitute}

\title{reconstitute}
\description{reconstitute}
\usage{
	reconstitute(d, conc=NULL, names=NULL)
	is.TSdatareconstitute(x)
}
\arguments {
	\item{d}{A concentrated data object.}
	\item{conc}{A concentrator.}
	\item{names}{Series names for the result.}
}
\value{Depends on the argument.}
\details{
	A concentrated data object is used to reconstruct the full dimension 
	data. Thus the result has the same dimension as the original data, but
	will not be exactly teh same because some information is lost when the
	data is concentrated. 
}
\seealso{
	\code{\link{concentrate}}
}
\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	z <- reconstitute(z)
	is.TSdatareconstitute(z)
}
\keyword{DSE}
\keyword{time series}


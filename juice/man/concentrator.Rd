\name{concentrator}
\alias{concentrator}
\alias{concentrator.concentrate}
\alias{concentrator.concentrator}
\alias{concentrator.TSdata}
\alias{concentrator.TSdataconcentrate}
\alias{concentrator.TSdataconcentrator}
\alias{concentrator.TSmodelconcentrate}
\alias{is.concentrator}
\alias{is.TSdataconcentrator}

\title{Data Concentrator}
\description{concentrator}
\usage{
	concentrator(d, ...)
	is.concentrator(x)
	is.TSdataconcentrator(x)
}
\arguments {
	\item{d}{A concentrate or concentrator object.}
}
\value{A concentrator.}
\details{
	The concentrator is extracted from a concentrated data object.
}
\seealso{
	\code{\link{concentrate}}
	\code{\link{concentrateOnly}}
	\code{\link{concentrateOriginal}}
}

\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	concentrator(z)
	is.concentrator(concentrator(z))
}
\keyword{DSE}
\keyword{time series}


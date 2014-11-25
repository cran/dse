\name{concentrateOriginal}
\alias{concentrateOriginal}
\alias{concentrateOriginal.TSdataconcentrate}
\alias{concentrateOriginal.TSdatareconstitute}
\alias{concentrateOriginal.concentrate}
\alias{concentrateOriginal.TScanonical.prediction}

\title{Extract Original Series}
\description{concentrateOriginal}
\usage{
	concentrateOriginal(d, ...)
}
\arguments {
	\item{d}{A concentrate object.}
}
\value{A TSdata object.}
\details{
	The original data set is returned as a TSdata object, stripped of
	the fact that it is a concentrate.
}
\seealso{
	\code{\link{concentrate}}
	\code{\link{concentrator}}
	\code{\link{concentrateOnly}}
}

\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	z <- concentrateOriginal(z)
}
\keyword{DSE}
\keyword{time series}


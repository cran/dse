\name{est.concentrated.model}
\alias{est.concentrated.model}
\alias{est.concentrated.model.TSdata}
\alias{est.concentrated.model.TSdataconcentrate}

\title{Estimate a Concentrated Model}
\description{est.concentrated.model}
\usage{
	est.concentrated.model(data, estimation="est.VARX.ls",
		estimation.args=NULL, m=1,p=1, center=TRUE, scale=TRUE)
	est.concentrated.model(data, estimation="est.VARX.ls",
		estimation.args=NULL)
}
\arguments {
	\item{data}{A TSdata or TSdataconcentrate object.}
	\item{estimation}{Estimation method.}
	\item{estimation.args}{Estimation method arguments.}
}
\value{A TSmodelconcentrate.}
\details{
	A concentrated version of the data (reduced dimension) is used to 
	estimate a reduced dimesion model. The projections for concentrating the
	data are retained so that model predictions can be expanded to the full
	dimesion data space.
	
	If data is TSdataconcentrate then the concentrator with that data is 
	used and m, p, center and scale are not used. For TSdata these arguments
	are used to first estimate a concentrated version of the data.
}
\seealso{
	\code{\link{concentrate}}
	\code{\link{concentrator}}
	\code{\link{est.projection}}
}

\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.concentrated.model(eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


\name{concentrate}
\alias{concentrate}
\alias{est.concentrate}
\alias{concentrate.default}
\alias{concentrate.TSdata}
\alias{is.concentrate}
\alias{is.TSdataconcentrate}
\alias{is.TSmodelconcentrate}

\title{Concentrate Series in a TSdata Object}
\description{
	Calculate a reduced dimesion version of the data using principal 
	components (or cannonical correlation for TSdata with input and output).
}
\usage{
   concentrate(data, conc=est.projection(data, center=center, scale=scale, n=n))
   concentrate(data, center=TRUE, scale=TRUE, n=1)     # data a matrix
   concentrate(data, center=TRUE, scale=TRUE, m=1,p=1) # data a TSdata object
   is.concentrate(x)
   is.TSdataconcentrate(x)
   is.TSmodelconcentrate(x)
   }
\arguments{
	\item{data}{A matrix or TSdata object.}
	\item{conc}  {object containing the projection matrix used for the reduction}
	\item{center}{center the observations to mean zero first (passed to est.projection).}
	\item{scale} {scale the observations to SD one first (passed to est.projection).}
	\item{n, m, p}{dimension of the concentrated series (passed to est.projection).}
}
\value{A matrix or TSdata object.}
\seealso{
	\code{\link{est.projection}}
	\code{\link{prcomp}}
}
\examples{
	if(is.R()) {
	   data("eg1.DSE.data", package="dse1")
	   require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	is.concentrate(z)
}
\keyword{DSE}
\keyword{time series}


\name{est.projection}
\alias{est.projection}
\alias{est.projection.default}
\alias{est.projection.TSdata}
\title{Calculate Projection from Concentrating Series.}
\description{
Calculate the projection to use for a reduced dimesion version of the data 
using principal components (or
cannonical correlation for TSdata with input and output).
}
\usage{
   est.projection(data, center=TRUE, scale=TRUE, n=1)     # data a matrix
   est.projection(data, center=TRUE, scale=TRUE, m=1,p=1) # data a TSdata object
   }
\arguments{
\item{data}{A matrix or TSdata object.}
\item{n, m, p}{dimension of the concentrated series.}
\item{center}{center the observations to mean zero first.}
\item{scale} {scale the observations to SD one first.}
}
\value{An object containing matrix (conc) to use to 
 concentrate the data.}
\seealso{
\code{\link{concentrate}}
\code{\link{reconstitute}}
\code{\link{prcomp}}
}
\examples{
if(is.R()) {
   data("eg1.DSE.data", package="dse1")
   require("mva")
  }
z <- est.projection(eg1.DSE.data)
}
\keyword{DSE}
\keyword{time series}


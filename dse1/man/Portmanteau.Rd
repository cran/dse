\name{Portmanteau}
\alias{Portmanteau}
\title{Calculate Portmanteau statistic}
\description{Calculate Portmanteau statistic.}
\usage{
	Portmanteau(res)}
\arguments{
\item{model}{A matrix with time-series residuals in columns.}
}
\seealso{
\code{\link{information.tests}}
}
\examples{
if(is.R()) require("ts")
Portmanteau(matrix(rnorm(200), 100,2)) # but typically with a residual
}
\keyword{DSE}
\keyword{ts}


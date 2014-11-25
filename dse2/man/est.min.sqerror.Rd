\name{est.min.sqerror}
\alias{est.min.sqerror}
\title{Estimate Min Squared Error}
\description{est.min.sqerror}
\usage{
	est.min.sqerror <- function(data, init.model, error.weights=1, ...)
}
\arguments {
	\item{data}{A TSdata object.}
	\item{init.model}{A starting model for minimization.}	
	\item{error.weights}{Weight to place on periods in the forecast.}
}
\value{A TSestModel}
\details{
	The function minimize is called with sum.sqerror as the objective
	function and error.weights as the argument. error.weights indicates 
	the weight to place on periods in the forecast.
 
}
\seealso{
	\code{\link{minimize}}
}

\keyword{DSE}
\keyword{ts}


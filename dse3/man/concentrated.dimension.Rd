\name{concentrated.dimension}
\alias{concentrated.dimension}
\alias{concentrated.dimension.concentrate}

\title{Concentrated Dimension}
\description{
	The dimension (number of series) in concentrate data. This is the
	dimension onto which the original series has been projected.
}
\usage{
	concentrated.dimension(x)
	concentrated.input.dimension(x)
	concentrated.output.dimension(x)
}
\arguments {
	\item{x}{A concentrated data object.}
}
\value{
	Depends on the argument. For a simple concentrated data object 
	the result is an integer.
}

\seealso{
	\code{\link{concentrated.input.dimension}}
	\code{\link{concentrated.output.dimension}}
	\code{\link{concentrate}}
}
\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	concentrated.dimension(output.data(z))
	concentrated.output.dimension(z)
}
\keyword{DSE}
\keyword{time series}


\name{concentrated.input.dimension}
\alias{concentrated.input.dimension}
\alias{concentrated.output.dimension}
\title{Concentrated Dimension of TSdata}
\description{
	The dimension (number of series) in concentrate data. This is the
	dimension onto which the original series has been projected.
}
\usage{
	concentrated.input.dimension(x)
	concentrated.output.dimension(x)
}
\arguments {
	\item{x}{A concentrated TSdata object.}
}
\value{An integer.}

\seealso{
	\code{\link{concentrated.dimension}}
	\code{\link{concentrate}}
}
\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	concentrated.output.dimension(z)
}
\keyword{DSE}
\keyword{time series}


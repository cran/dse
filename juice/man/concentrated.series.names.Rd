\name{concentrated.series.names}
\alias{concentrated.series.names}
\alias{concentrated.series.names.concentrate}
\alias{concentrated.series.names.TSdata}
\alias{concentrated.input.series.names}
\alias{concentrated.output.series.names}

\title{Concentrated Series Names}
\description{
	The names of series in concentrate data.
}
\usage{
	concentrated.series.names(x)
	concentrated.input.series.names(x)
	concentrated.output.series.names(x)
}
\arguments {
	\item{x}{A concentrated data object.}
}
\value{
	Depends on the argument. For a simple concentrated data object 
	the result is a vector of strings.
}

\seealso{
	\code{\link{series.names}}
	\code{\link{input.series.names}}
	\code{\link{output.series.names}}
	\code{\link{concentrated.dimension}}
	\code{\link{concentrate}}
}
\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	concentrated.series.names(z)
}
\keyword{DSE}
\keyword{time series}


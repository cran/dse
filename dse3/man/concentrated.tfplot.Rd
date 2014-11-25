\name{concentrated.tfplot}
\alias{concentrated.tfplot}
\title{Plot Concentrated Series}
\description{
	The concentrate data is plotted.
}
\usage{
	concentrated.tfplot(x)
}
\arguments {
	\item{x}{A concentrated data object.}
}
\value{
	Depends on the argument. For a simple concentrated data object 
	the result is a vector of strings.
}
\sideffects{A plot is generated}

\seealso{
	\code{\link{tfplot}}
}
\examples{
	if(is.R()) {
	data("eg1.DSE.data", package="dse1")
	require("mva")
	}
	z <- concentrate(eg1.DSE.data)
	concentrated.tfplot(z)
}
\keyword{DSE}
\keyword{time series}


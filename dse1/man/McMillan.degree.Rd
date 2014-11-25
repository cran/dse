\name{McMillan.degree}
\alias{McMillan.degree}
\alias{McMillan.degree.ARMA}
\alias{McMillan.degree.SS}
\alias{McMillan.degree.TSestModel}

\title{calculate McMillan degree}
\description{

Calculate the McMillan degree of an ARMA TSmodel.
}
\usage{
	McMillan.degree(model, fuzz=fuzz, verbose=T, warn=T)
}
\arguments{
\item{ model}{An object of class TSmodel.}

\item{ fuzz}{Roots within fuzz distance are counted as equivalent.}
\item{ verbose}{If T roots are printed.}
\item{ warn}
{If F then warnings about unit roots added for TREND are not printed.}
}
\value{
A list with elements gross and distinct containing all roots 
and distinct roots.
}
\sideffects{
The number of roots and distinct roots is printed if verbose is T.
}
\seealso{

\code{\link{stability}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
McMillan.degree(model)
}
\keyword{DSE}
\keyword{ts}


\name{polynomials}
\alias{characteristic.poly}
\alias{companion.matrix}
\alias{poly.value}
\alias{polydet}
\alias{polyprod}
\alias{polysum}
\alias{polyroot.det}
\alias{old.polyprod}
        
\title{Polynomial Utilities}
\description{Polynomial utility functions used by DSE.}
\usage{
	characteristic.poly(x)
	companion.matrix(x)
	poly.value(x)
	polydet(x)
	polyprod(x,y)
	polysum(x,y)
	polyroot.det(x)
	old.polyprod(x,y)
}
\arguments {
	\item{x}{An array representing a polynomial.}
	\item{y}{An array representing a polynomial.}
}
\value{depends}
\details{
	These are utility functions used in some ARMA model calculations
	such as root and stability calculations.}
\seealso{
	\code{\link{polyroot}}
	\code{\link{roots}}
	\code{\link{stability}}
}

\keyword{DSE}
\keyword{ts}


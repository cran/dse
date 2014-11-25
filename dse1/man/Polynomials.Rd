\name{Polynomials}
\alias{characteristic.poly}
\alias{companion.matrix}
\alias{poly.value}
\alias{polydet}
\alias{polyprod}
\alias{polysum}
\alias{polyroot.det}
        
\title{Polynomial Utilities}
\description{Polynomial utility functions used by DSE.}
\usage{
    characteristic.poly(A)
    companion.matrix(A)
    poly.value(coef, z)
    polydet(A)
    polyprod(a, b)
    polysum(a, b)
    polyroot.det(A)
}
\arguments{
    \item{A}{An array representing a matrix polynomial.}
    \item{B}{An array representing a matrix polynomial.}
    \item{a}{An array representing a matrix polynomial.}
    \item{b}{An array representing a matrix polynomial.}
    \item{coef}{Coefficients of a polynomial.}
    \item{z}{Value at which the polynomial is to be evaluated.}
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

%\keyword{DSE}
\keyword{ts}


\name{richardson.grad}
\alias{richardson.grad}
\title{Gradient of a Function}
\description{Calculate the gradient of a function by Richardson improvement.}
\usage{
	richardson.grad(func,x, d=0.01, eps=1e-4, r=6, show.details=F)
}
\arguments{
\item{ func}
{A function with a real result.}
\item{ x}
{
A real or real vector argument to func, indicating the point at which the
gradient is to be calculated.}

\item{ d eps}
{
The epsilon interval value used for the first Richardson iteration is set 
to d*x or eps if x is 0.0.}
\item{ r}{The number of Richardson improvement iterations.}
\item{ show.details}{If T intermediate calculations are printed.}
}
\value{A real or real vector of the calculated gradient.}
\seealso{\code{\link{numerical.grad}}}
\examples{richardson.grad(sin, pi)}
\keyword{DSE}
\keyword{ts}


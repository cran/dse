\name{numerical.grad}
\alias{numerical.grad}
\title{Gradient of a Function}
\description{
Calculate the gradient of a function by simple epsilon differencing.
}
\usage{
	numerical.grad(func,x, eps=1e-12)
}
\arguments{
\item{func}{A function with a real result.}
\item{x}{
A real or real vector argument to func, indicating the point at which the
gradient is to be calculated.
}

\item{eps}
{The epsilon difference to use for calculating the gradient.}
}
\value{A real or real vector of the calculated gradient.}
\seealso{\code{\link{richardson.grad}}}
\examples{numerical.grad(sin, pi)}
\keyword{DSE}
\keyword{ts}


\name{span}
\alias{span}
\alias{span.v0}
\alias{span.v00}
\alias{span.default}
\alias{span.TSestModel}

\title{Calculate Span of Tangent Plan}
\description{Calculate the dimension of the tangent space}
\usage{
	span(func, x,d=0.01, eps=1e-4,r=6, show.details=F)
}
\arguments{
\item{func}{An S function. func must have a single vector argument.}
\item{x}{The parameter vector argument to func.}
\item{d}{The fraction of x to use for the initial numerical approximation.}
\item{eps}{Used instead of d for elements of x which are zero.}
\item{r}{The number of Richardson improvement iterations.}
}
\value{
The singular values of the matrix of tangent vectors are returned.
}
\sideffects{
If show.details is T then intermediate calculations are printed.
}
\details{
span performs a svd of the tangent vectors at the point x. This can be used
to calculate the dimension of the tangent space (ie. by over specifying 
the model and counting the number of significant singular values).
This function uses Richardson extrapolation  (for more details      
see the functions richardson.grad and genD) to get a numerical 
approximation of the tangent vectors to the parameter 
manifold. SVD is then used to calculate their span.
}
\seealso{
\code{\link{richardson.grad}},
\code{\link{genD}}
}
\examples{
func <- function(x){c(x[1], x[1], x[2]^2)}
span("func", c(2,2))
span("func", c(2,5))
span("func", c(2,2,5))
}
\keyword{DSE}
\keyword{ts}


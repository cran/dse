\name{genD}
\alias{genD}
\alias{genD.default}
\alias{genD.default.v2}
\alias{genD.ARMA}
\alias{genD.innov}
\alias{genD.TSestModel}
\alias{genD.c}

\title{genD}
\description{Generate a matrix of function derivative information.}
\usage{
	genD(func, x, d=0.01, eps=1e-4, r=6)
	genD.c(func, x, d=0.01, eps=1e-4, r=6)
	genD.TSestModel(estModel, d=0.01, eps=1e-4, r=6)
}
\arguments{
\item{ func, model}
{An S function or an object of class TSestModel. func must have a single
vector argument, in which case the second argument,x, is the parameter
vector.}
\item{d}{The fraction of x to use for the initial numerical approximation.}
\item{eps}{Used instead of d for elements of x which are zero.}
\item{r}{The number of Richardson improvement iterations.}
}
\value{
A list with three elements as follows:
   $D is a matrix of first(gradients) and second order partial
      derivatives organized in the same manner as Bates and 
      Watts. (The first p columns are the gradients and the 
      next p(p-1)/2 columns are the lower triangle of the Hessian).
   $p is the dimension of the parameter space=dim of the tangent space.
   $f0 is the function value at the point where the matrix D 
        was calculated. 
}
\details{The derivatives are calculated numerically using Richardson improvement.}
\seealso{
\code{\link{curvature}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(TSdata( output=output.data(eg1.DSE.data.diff, series=1:2)),
                      max.lag=2)
z <- genD(model)
}
\keyword{DSE}
\keyword{ts}


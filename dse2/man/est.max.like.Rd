\name{est.max.like}
\alias{est.max.like}
\title{Maximum Likelihood Estimation}
\description{Maximum likelihood estimation.}
\usage{est.max.like(model,  max.iter=20, ftol=1e-5, algorithm='nlm')}
\arguments{
\item{model}
{An object of class TSestModel, or a previously returned result from max.like.}
\item{ max.iter}{The maximum number of iterations.}
\item{ ftol}{The function tolerance for indicating convergence.}
\item{ algorithm}{The algorithm ('nlmin' or 'dfpMin') to use for maximization.}
}
\value{
The value returned is an object of class TSestModel with  additional
elements \code{converged}, which is T or F indicating convergence, 
and \code{dfpMin.results} or \code{nlmin.results}.
If this function calls dfp the Hessian,etc are return as \code{dfpMin.results}.
If this function is called again and those results are
available then they are used. 
This could cause problems if element \code{model} is modified. If that is
done then the element \code{dfpMin.results} should be set to NULL.
}

\note{THE CURRENT VERSION OF THIS PROGRAM IS SUPERFICIALLY BROKEN.}

\seealso{
\code{\link{nlm}}
\code{\link{est.VARX.ls}}
\code{\link{bft}}
}
\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
# est.max.like(est.VARX.ls(egJofF.1dec93.data)) DSEbug needs work
}
\keyword{DSE}
\keyword{ts}


\name{distribution}
\alias{distribution}
\alias{distribution.default}
\alias{distribution.TSdata}
\alias{distribution.estimation.evaluation}
\alias{distribution.parms.ee}
\alias{distribution.roots.ee}
\title{Plot distribution of estimates}
\description{Plot distribution of estimates.}
\usage{
	distribution(obj,  Sort=F, bandwidth=0.2)
}
\arguments{
\item{ obj}{An object as returned by eval.estimation.} 
\item{ Sort}
{If Sort is true then sort is applied. This helps (a bit)
with estimation methods like black.box which may not return parameters
of the same length or in the same order.}
\item{ bandwidth}{Passed to \code{density} or \code{ksmooth}.}
}
\value{None}
\details{
\code{ksmooth} is applied if available to get a smoothed 
estimate of the distribution of the estimates. If \code{ksmooth} is not
available then \code{density} is applied if it is available.
}
\seealso{
\code{\link{eval.estimation}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(TSdata(output=output.data(eg1.DSE.data.diff)))
z <-  eval.estimation(model, 
	estimation="est.VARX.ls", estimation.args=list(max.lag=2))
#distribution(z) Rbug v90.1
}
\keyword{DSE}
\keyword{ts}


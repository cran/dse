\name{reachability}
\alias{reachability}
\alias{reachability.ARMA}
\alias{reachability.SS}
\alias{reachability.TSestModel}

\title{Calculate reachability of a TSmodel}
\description{
	Calculate the singular values of the reachability matrix of a model.
}
\usage{
	reachability(model)
}
\arguments{
\item{ model}{An object containing TSmodel.}}
\value{The singular values of the reachability matrix.}
\details{

If all singular values are significantly different from zero the
model is controllable.
}
\seealso{

\code{\link{observability}},
\code{\link{stability}}
\code{\link{roots}}
\code{\link{McMillan.degree}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
reachability(model)
}
\keyword{DSE}
\keyword{ts}


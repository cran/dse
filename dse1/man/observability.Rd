\name{observability}
\alias{observatility}
\alias{observability.ARMA}
\alias{observability.SS}
\alias{observability.TSestModel}

\title{Calculate observability of a TSmodel}
\description{
	Calculate the singular values of the observability matrix of a model.
}
\usage{
	observability(model)
}
\arguments{
	\item{model}{An object containing a TSmodel.}
}
\value{The singular values of the observability matrix.}
\details{
	If all singular values are significantly different from zero the
	model is observable.
}
\seealso{
	\code{\link{reachability}},
	\code{\link{stability}}
	\code{\link{McMillan.degree}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
	observability(model)
}
\keyword{DSE}
\keyword{ts}


\name{to.SS.innov}
\alias{to.SS.innov}
\title{Convert to State Space Innovations Model}
\description{Convert to a state space innovations representation.}
\usage{
	to.SS.innov(model)}
\arguments{
	\item{model}{An object of class TSmodel.}
}
\value{
An object of class 'SS' 'TSmodel' containing an 
innovations state space representation.
}
\seealso{
\code{\link{to.SS}},
\code{\link{to.SS.Oform}}
\code{\link{to.SS.Chol}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
model <- to.SS.innov(model)
}
\keyword{DSE}
\keyword{ts}


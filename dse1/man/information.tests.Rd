\name{information.tests}
\alias{information.tests}
\title{Tabulates selection criteria}
\description{Tabulates several model selection criteria.}
\usage{
	information.tests(..., sample.start=1,sample.end=NULL, Print=T, warn=T)}
\arguments{
\item{...}{At least one object of class TSestModel.}
\item{sample.start}{The start of the period to use for criteria calculations.}
\item{sample.end}  {
	The end of the period to use for criteria calculations. If omitted
	the end of the sample is used.}
}
\item{Print}{If FALSE then printing suppressed.}
\item{warn} {If FALSE then some warning messages are suppressed.}

\value{A matrix of the value for each model on each test returned invisibly.}
\sideffects{Criteria are tabulated for all models in the list.}
\seealso{
	\code{\link{information.tests.calculations}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model1 <- est.VARX.ls(eg1.DSE.data.diff)
	model2 <- est.VARX.ar(eg1.DSE.data.diff)
	information.tests(model1, model2)}
\keyword{DSE}
\keyword{ts}


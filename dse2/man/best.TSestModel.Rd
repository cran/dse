\name{best.TSestModel}
\alias{best.TSestModel}
\title{Select Best Model}
\description{Select the best model.}
\usage{
	best.TSestModel(models, sample.start=10, sample.end=NULL,
	 criterion='aic', verbose=T)
}
\arguments{
\item{ models}{A list of TSestModels.}
\item{ sample.start}{The starting point to use for 
calculating information criteria.}
\item{ criteria}
{The selection criteria. 'taic' would be a better choice but this is not
available for VAR and ARMA models.}
}
\value{A TSestModel}
\details{Information criteria are calculated and the best model returned.}
\seealso{
\code{\link{est.black.box1}},
\code{\link{est.black.box2}}
\code{\link{est.black.box3}}
\code{\link{est.black.box4}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
models <- list(est.VARX.ls(eg1.DSE.data.diff), est.VARX.ar(eg1.DSE.data.diff))
z <-  best.TSestModel(models)
}
\keyword{DSE}
\keyword{ts}


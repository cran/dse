\name{best.TSestModel}
\alias{best.TSestModel}
\title{Select Best Model}
\description{Select the best model.}
\usage{
    best.TSestModel(models, sample.start=10, sample.end=NULL,
     criterion='aic', verbose=TRUE)
}
\arguments{
    \item{models}{a list of TSestModels.}
    \item{sample.start}{the starting point to use for  calculating 
       information criteria.}
    \item{sample.end}{the end point to use for  calculating 
       information criteria.}
    \item{criterion}{Criterion to be used for model 
       selection. see \code{information.tests.calculations}. 'taic' would 
       be a better default
       but this is not available for VAR and ARMA models.}
    \item{verbose}{if TRUE then additional information is printed.}
}
\value{A TSestModel}
\details{Information criteria are calculated and the best model returned.}
\seealso{
\code{\link{est.black.box1}},
\code{\link{est.black.box2}}
\code{\link{est.black.box3}}
\code{\link{est.black.box4}}
\code{\link[dse1]{information.tests.calculations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
models <- list(est.VARX.ls(eg1.DSE.data.diff), est.VARX.ar(eg1.DSE.data.diff))
z <-  best.TSestModel(models)
}
%\keyword{DSE}
\keyword{ts}


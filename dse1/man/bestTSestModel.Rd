\name{bestTSestModel}
\alias{bestTSestModel}
\title{Select Best Model}
\description{Select the best model.}
\usage{
    bestTSestModel(models, sample.start=10, sample.end=NULL,
     criterion='aic', verbose=TRUE)
}
\arguments{
    \item{models}{a list of TSestModels.}
    \item{sample.start}{the starting point to use for  calculating 
       information criteria.}
    \item{sample.end}{the end point to use for  calculating 
       information criteria.}
    \item{criterion}{Criterion to be used for model 
       selection. see \code{informationTestsCalculations}. 'taic' would 
       be a better default
       but this is not available for VAR and ARMA models.}
    \item{verbose}{if TRUE then additional information is printed.}
}
\value{A TSestModel}
\details{Information criteria are calculated and the best model returned.}
\seealso{
\code{\link{estBlackBox1}},
\code{\link{estBlackBox2}}
\code{\link{estBlackBox3}}
\code{\link{estBlackBox4}}
\code{\link[dse1]{informationTestsCalculations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
models <- list(estVARXls(eg1.DSE.data.diff), estVARXar(eg1.DSE.data.diff))
z <-  bestTSestModel(models)
}
%\keyword{DSE}
\keyword{ts}


\name{est.black.box1}
\alias{est.black.box1}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    est.black.box1(data, estimation="est.VARX.ls", 
        reduction="reduction.Mittnik", 
        criterion="taic", trend=F, subtract.means=F, verbose=T, max.lag=6)
}
\arguments{
    \item{data}{Data in an object of class TSdata.}
    \item{estimation}{Initial estimation method to be used.}
    \item{reduction}{Reduction method to be used.}
    \item{criterion}{
    Criterion to be used for model selection. This can be: 
port, like, aic, bic, gvc, rice,f pe, taic, tbic, tgvc, trice 
or tfpe for the  Portmanteau test, likelihood, Akaike Information Criterion, 
Bayes  Information Criterion, Generalized Cross Validation, Rice Criterion, 
or Final Prediction Error. The preceeding 't' indicates that the theoretical 
parameter space dimension should be used. }
    \item{verbose}{Indicates if information should be printed during estimation.}
}
\value{A state space model in an object of class TSestModel.}
\sideffects{
If verbose is T then estimation information is printed and check.residuals is run, 
which gives plots of information about the residuals.
}
\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
goodmodel <- est.black.box1(egJofF.1dec93.data)
}
%\keyword{DSE}
\keyword{ts}


\name{estBlackBox1}
\alias{estBlackBox1}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    estBlackBox1(data, estimation="estVARXls", 
        reduction="MittnikReduction", 
        criterion="taic", trend=FALSE, subtract.means=FALSE, 
	verbose=TRUE, max.lag=6)
}
\arguments{
    \item{data}{Data in an object of class TSdata.}
    \item{estimation}{Initial estimation method to be used.}
    \item{reduction}{Reduction method to be used.}
    \item{criterion}{Criterion to be used for model 
       selection. see \code{informationTestsCalculations}.}
    \item{trend}{logical indicating if a trend should be estimated.}
    \item{subtract.means}{logical indicating if the mean should be subtracted
       from data before estimation.}
    \item{verbose}{logical indicating if information should be printed 
       during estimation.}
    \item{max.lag}{integer indicating the maximum number of lags to consider.}
}
\value{A state space model in an object of class TSestModel.}
\sideffects{
If verbose is TRUE then estimation information is printed and checkResiduals is run, 
which gives plots of information about the residuals.
}
\seealso{
    \code{\link[dse1]{informationTestsCalculations}}
}

\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
goodmodel <- estBlackBox1(egJofF.1dec93.data)
}
\concept{DSE}
\keyword{ts}


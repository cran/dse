\name{est.VARX.ar}
\alias{est.VARX.ar}
\alias{old.est.VARX.ar}

\title{Estimate a VAR TSmodel}
\description{Estimate a VAR TSmodel with (optionally) an exogenous input.}
\usage{
    est.VARX.ar(data, subtract.means=FALSE,  re.add.means=TRUE, standardize=FALSE, 
         unstandardize=TRUE, aic=TRUE, max.lag=NULL, method="yule-walker", warn=TRUE)
}
\arguments{
    \item{data}{A TSdata object.}
    \item{subtract.means}{
        If TRUE subtract the means from the data before estimation.}
    \item{re.add.means}{
      If TRUE the model is adjusted for the non-zero mean data when returned. 
      If subtract.means is also TRUE then the mean is added back to the data.}
    \item{standardize}{
      Note that the mean is not subtracted unless subtract.means is TRUE.
      A VAR model in an object of class TSestModel.}
    \item{unstandardize}{
      If TRUE and standardize is TRUE then the returned model is adjusted to 
      correspond to the original data.}
    \item{aic}{Passed to function ar.}
    \item{max.lag}{The maximum number of lags that should be considered.}
    \item{method}{Passed to function ar.}
    \item{warn}{If TRUE certain warning message are suppressed.}
}
\value{A TSestModel object containing a TSmodel object which is aVAR model.}

\details{
This function uses the function ar. The ar 
method does not support trend estimation (as in est.VARX.ls).
If aic=TRUE the number of lags is determined by an AIC statistic (see ar). 
If an exogenous (input) 
variable is supplied the input and output are combined (i.e.- both
treated as outputs) for estimation, and the resulting model is
converted back by transposing the exogenous variable part of the 
polynomial and discarding inappropriate blocks.
Residuals,etc, are calculated by evaluating the estimated model as a 
TSmodel/ARMA with the data (ie. residuals are not the residuals from the
regression).
Note: ar uses a Yule-Walker approach (uses autocorrelations) so effectively the 
  model is for data with means removed. Thus subtract.means does not make much
  difference and re.add.means must be TRUE to get back to a model for the 
  original data.
}
\references{
    Gilbert, P. D. (1993) State space and ARMA models: An overview of
    the equivalence. Working paper 93-4, Bank of Canada. Available at
    <www.bank-banque-canada.ca/pgilbert>

    Gilbert, P. D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.
}
\seealso{
    \code{\link{est.SS.from.VARX}}
    \code{\link[dse2]{est.SS.Mittnik}}
    \code{\link[dse2]{bft}}
    \code{\link{est.VARX.ls}}
    \code{\link[ts]{ar}}
    \code{\link{DSE.ar}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ar(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


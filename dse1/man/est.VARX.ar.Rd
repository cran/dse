\name{est.VARX.ar}
\alias{est.VARX.ar}
\alias{old.est.VARX.ar}

\title{Estimate a VAR TSmodel}
\description{Estimate a VAR TSmodel with (optionally) an exogenous input.}
\usage{
	est.VARX.ar(data,subtract.means=F,re.add.means=T, 
		standardize=F, unstandardize=T, aic=T, max.lag=NULL, 
		method='yule-walker')
}
\arguments{
	\item{data}{A TSdata object.}
	\item{subtract.means}
	  {If T subtract the means from the data before estimation.}
	\item{re.add.means}{
	  If T the model is adjusted for the non-zero mean data when returned. 
	  If subtract.means is also T then the mean is added back to the data.}
	\item{standardize}{
	  Note that the mean is not subtracted unless subtract.means is T.
	  A VAR model in an object of class TSestModel.}
	\item{unstandardize}{
	  If T and standardize is T then the returned model is adjusted to 
	  correspond to the original data.}
}
\value{A TSestModel object containing a TSmodel object which is aVAR model.}

\details{
This function uses the Splus function ar. The ar 
method does not support trend estimation (as in est.VARX.ls).
If aic=T the number of lags is determined by an AIC statistic (see ar). 
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
  difference and re.add.means must be T to get back to a model for the 
  original data.
}
\section{ REFERENCES} {
	P.D.Gilbert (1993) State Space and ARMA models: An overview of
	estimation and reduction.
}
\seealso{
	\code{\link{est.SS.from.VARX}}
	\code{\link[dse2]{est.SS.Mittnik}}
	\code{\link[dse2]{bft}}
	\code{\link{est.VARX.ls}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ar(eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


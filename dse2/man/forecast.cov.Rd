\name{forecast.cov}
\alias{forecast.cov}
\alias{forecast.cov.TSdata}
\alias{forecast.cov.TSmodel}
\alias{forecast.cov.TSestModel}
\alias{forecast.cov.compiled}
\alias{forecast.cov.compiled.ARMA}
\alias{forecast.cov.compiled.SS}
\alias{forecast.cov.compiled.innov}
\alias{forecast.cov.compiled.non.innov}
\alias{forecast.cov.reductions.wrt.true}
\alias{forecast.cov.single.TSmodel}
\alias{is.forecast.cov}

\title{Forecast covariance for different models}
\description{Calculate the forecast covariance for different models.}
\usage{
	forecast.cov(obj, ..., data=NULL, discard.before=NULL,
	   horizons=1:12, zero=F, trend=F, estimation.sample= periods(data))
	is.forecast.cov(obj)
}
\arguments{
	\item{obj}{TSdata or one or more TSmodels or TSestModels}
	\item{data}{An object of class TSdata.}
	\item{discard.before}{
	Period before which forecasts should be discarded when calculating covariance.}
	\item{horizons}{Horizons for which forecast covariance should be calculated.}
	\item{zero}{If T the covariance is calculated for a forecast of zero.}
	\item{trend}{If T the covariance is calculated for a forecast of trend.}
	\item{estimation.sample}{Portion of the sample to use for calculating the trend.}
}
\value{
	A list with the forecast covariance for supplied models on the given
	sample. This is in the element \code{forecast.cov} of the result. Other
	elements contain information in the arguments.
}
\details{
Calculate the forecast cov of obj relative to data. If obj is TSdata
then the output data is used as the forecast. For other classes of obj
TSmodel(obj) is used with data to produce a forecast. TSmodel() is also
applied to each element of ... to extract a model. All models should
work with data. If obj is a TSestModel and data is NULL then TSdata(obj)
is used as the data.
This is just multiple applications of forecast.cov.single.TSmodel
discard.before is an integer indicating the number of points in the
beginning of forecasts to discard before calculating covariances.
If it is the default, NULL, then the default (minimum.startup.lag) will
be used for each model and the default (1) will be used for trend and zero.
If zero  is T then forecast.cov is also calculated for a forecast of zero.
If trend is T then forecast.cov is also calculated for a forecast of a linear 
trend using data to estimation.sample.
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model1 <- est.VARX.ar(eg1.DSE.data.diff)
	model2 <- est.VARX.ls(eg1.DSE.data.diff)
	z <-  forecast.cov(model1, model2, data=trim.na(eg1.DSE.data.diff))
	is.forecast.cov(z)
}
\keyword{DSE}
\keyword{ts}


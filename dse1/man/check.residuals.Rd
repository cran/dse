\name{check.residuals}
\alias{check.residuals}
\alias{check.residuals.default}
\alias{check.residuals.TSdata}
\alias{check.residuals.TSestModel}

\title{Autocorrelations Diagnostics}
\description{Calculate autocorrelation diagnostics of a time series 
matrix or TSdata or residuals of a TSestModel}
\usage{
	check.residuals(model, ac=T, pac=T, select=NULL, drop=NULL, 
	  plot.=T,verbose=F)
}
\arguments{
	\item{model}{An TSestModel or TSdata object.}
	\item{ac}{If T the auto-correlation function is plotted.}
	\item{pac}{If T the partial auto-correlation function is plotted.}
	\item{select}{
	  Is used to indicate a subset of the residual series. By default all
	  residuals are used.}
	\item{drop}{
	  Is used to indicate a subset of the residual time periods to drop. 
	  All residuals are used with the default (NULL).Typically this can be used 
	  to get rid of bad initial conditions (eg. drop=seq(10) ) or outliers.}
	\item{plot.}{If F then plots are not produced.}
	\item{verbose}{
		If T then the auto-correlations and partial auto-correlations 
		are printed if they are calculated.}
}
\value{
	A list with residual diagnostic information: residuals, mean, cov, 
	acf= autocorrelations, pacf= partial autocorrelations.
}
\details{
This is a generic function. The default method works for a time series
matrix which is treated as if it were a matrix of residuals. However, in a
Box-Jenkins type of analysis the matrix may be data which is being evaluated to
determine a model. The method for a TSestModel evaluates the residuals
calculated by subtracting the output data from the model predictions.}

\sideffects{
	Diagnostic information is printed and plotted if a device is available.
	Use dev.ask(T) to pause plotting between pages.
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.VARX.ls(eg1.DSE.data.diff)
	check.residuals(model)
}
\keyword{DSE}
\keyword{ts}


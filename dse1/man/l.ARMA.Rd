\name{l.ARMA}
\alias{l.ARMA}
\title{Evaluate an ARMA TSmodel}
\description{Evaluate an ARMA TSmodel.}
\usage{
l.ARMA(model, data, sampleT=NULL, predictT=NULL,result=NULL, compiled=.DSECOMPILED, warn=T)
}
\arguments{
\item{model}{An 'ARMA' 'TSmodel' object.}
\item{data}{A TSdata object.}
\item{sampleT} {An integer indicating the number of periods of data to use.}
\item{predictT}{An integer to what period forecasts should be extrapolated.}
\item{result}
{ If non-NULL then the returned value is only the sub-element indicated by 
result. result can be a character string or integer.}
\item{compiled}
{
Indicates if a call should be made to the compiled code for computation. 
A FALSE value is mainly for testing purposes.}
\item{warn}{If FALSE then certain warning messages are turned off.}
}
\value{An object of class TSestModel (see TSestModel).}

\details{
This function is called by the function l() when the argument to l is an ARMA
model. Using l() is usually preferable to calling l.ARMA directly.
l.ARMA calls a compiled program unless compiled=F. The compiled version is much 
faster. 
}
\seealso{
\code{\link{l}},
\code{\link{l.SS}}
\code{\link{TSmodel}}
\code{\link{TSestModel.object}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- TSmodel(est.VARX.ls(eg1.DSE.data.diff))
	evalutated.model <- l.ARMA(model,eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


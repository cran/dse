\name{l.SS}
\alias{l.SS}
\title{Evaluate a state space TSmodel}
\description{Evaluate a state space TSmodel.
}
\usage{
	l.SS( model, data, sampleT=NULL, predictT=NULL,
	   return.state=F, return.track=F,result=NULL, compiled=.DSECOMPILED)
}
\arguments{
\item{model}{An 'SS' 'TSmodel' object.}
\item{data}{A TSdata object.}
\item{sampleT}
{
An integer indicating the last data point to use for one step ahead filter
estimation. If NULL all available data is used. }
\item{predictT}
{
An integer indicating how far past the end of the sample predictions should 
be made. For models with an input, input data must be provided up to predictT. 
Output data is necessary only to sampleT. If NULL predictT is set to sampleT.}
\item{return.state}
{
If T the element $filter$state containing E[z(t)|y(t-1), u(t)] is returned as
part of the result. This can be a fairly large matrix. }
\item{return.track}
{If T the element $filter$track containing the expectation of the tracking 
 error given y(t-1) and u(t) is returned as part of the result. This can be an 
very large array. }
\item{result}
{If result is not specified an object of class TSestModel is 
returned. Otherwise, the specified element of TSestModel$estimates
is returned.}
\item{compiled}
{If T the compiled version of the code is used. Otherwise the S version is used.}
\item{warn}
{If FALSE then certain warning messages are turned off.}
}
\value{
Usually an object of class TSestModel (see TSestModel), but see result above.}

\details{
This function is called by the function l() when the argument to l is a state
space model. Using l() is usually preferable to calling l.SS directly.
l.SS calls a compiled program unless compiled=F. The compiled version is much faster
than the S version.
}
\seealso{
\code{\link{l}}
\code{\link{l.ARMA}}
\code{\link{TSmodel}}
\code{\link{TSestModel.object}}
\code{\link{smoother}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(TSmodel(est.VARX.ls(eg1.DSE.data.diff)))
lmodel <- l.SS(model,eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


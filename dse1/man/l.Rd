\name{l}
\alias{l}
\alias{l.TSdata}
\alias{l.TSestModel}

\title{Evaluate a TSmodel}
\description{Evaluate a model with data.}
\usage{
	l(model, data, sampleT=NULL, predictT=NULL,result=NULL)}
\arguments{
\item{model}{A TSmodel object.}
\item{data}{A TSdata object.}
\item{sampleT}
{
The length of data which should be used for estimation of one step ahead
predictions. Output data
must be at least as long as sampleT. If sampleT is not supplied it
is taken to be periods(data).}
\item{predictT}
{
The length of data which should be used for simulation (i.e. - multi-step 
prediction) beyond sampleT. Input data
must be at least as long as predictT. predictT must be at least as
large as sampleT. If predictT is not supplied it
is taken to be sampleT.}
\item{result}
{
If result is not specified an object of class TSestModel is 
returned. Otherwise, the specified element of the estimates
is returned.}
}
\value{See result above.}
\details{
This function determines whether the model is in state space or 
ARMA form and calls l.SS or l.ARMA.
}
\seealso{
\code{\link{l.SS}},
\code{\link{l.ARMA}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(TSmodel(est.VARX.ls(eg1.DSE.data.diff)))
evalutated.model <- l(model, eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


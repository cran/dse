\name{l.SS}
\alias{l.SS}
\title{Evaluate a state space TSmodel}
\description{Evaluate a state space TSmodel.
}
\usage{
    l.SS(obj1, data, sampleT=NULL, predictT=NULL, error.weights=0,
         return.state=FALSE, return.track=FALSE, result=NULL, 
	 compiled=.DSECOMPILED,
         warn=TRUE, return.debug.info=FALSE)
}
\arguments{
    \item{obj1}{An 'SS' 'TSmodel' object.}
    \item{data}{A TSdata object.}
    \item{sampleT}{an integer indicating the last data point to use 
       for one step ahead filter estimation. If NULL all available data is used. }
    \item{predictT}{an integer indicating how far past the end of the 
       sample predictions should be made. For models with an input, 
       input data must be provided up to predictT. Output data is 
       necessary only to sampleT. If NULL predictT is set to sampleT.}
    \item{error.weights}{a vector of weights to be applied to the 
       squared prediction errors.}
    \item{return.state}{if T the element \code{filter$state} 
       containing E[z(t)|y(t-1), u(t)] is returned as part of the 
       result. This can be a fairly large matrix.}
    \item{return.track}{if T the element \code{filter$track} containing
       the expectation of the tracking error given y(t-1) and u(t) is 
       returned as part of the result. This can be an very large array. }
    \item{result}{if result is not specified an object of class 
       TSestModel is returned. Otherwise, the specified element 
       of \code{TSestModel$estimates} is returned.}
    \item{compiled}{if T the compiled version of the code is used. 
       Otherwise the S/R version is used.}
    \item{warn}{if FALSE then certain warning messages are turned off.}
    \item{return.debug.info}{logical indicating if additional debugging
       information should be returned.}
}
\value{
Usually an object of class TSestModel (see TSestModel), but see result above.}

\details{
This function is called by the function l() when the argument to l is a state
space model. Using l() is usually preferable to calling l.SS directly.
l.SS calls a compiled program unless compiled=F. The compiled version is much faster
than the S version.

Output data must be at least as long as sampleT. If sampleT is not supplied it
is taken to be periods(data).

Input data must be at least as long as predictT. predictT must be at least as
large as sampleT. If predictT is not supplied it
is taken to be sampleT.

If \code{error.weights} is greater than zero then weighted prediction 
errors are calculated up to the horizon indicated
by the length of error.weights. The weights are applied to the squared
error at each period ahead.

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
%\keyword{DSE}
\keyword{ts}


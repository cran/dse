\name{l.ARMA}
\alias{l.ARMA}
\title{Evaluate an ARMA TSmodel}
\description{Evaluate an ARMA TSmodel.}
\usage{
l.ARMA(obj1, dat, sampleT=NULL, predictT=NULL,result=NULL,
       error.weights=0,  compiled=.DSECOMPILED, warn=TRUE, return.debug.info=FALSE)
}
\arguments{
    \item{obj1}{an 'ARMA' 'TSmodel' object.}
    \item{dat}{a TSdata object.}
    \item{sampleT}{an integer indicating the number of periods of data to use.}
    \item{predictT}{an integer to what period forecasts should be extrapolated.}
    \item{result}{
      if non-NULL then the returned value is only the sub-element indicated by 
      result. result can be a character string or integer.}
    \item{error.weights}{a vector of weights to be applied to the 
       squared prediction errors.}
    \item{compiled}{indicates if a call should be made to the compiled 
      code for computation. A FALSE value is mainly for testing purposes.}
    \item{warn}{if FALSE then certain warning messages are turned off.}
    \item{return.debug.info}{logical indicating if additional debugging
       information should be returned.}
}

\value{an object of class TSestModel (see TSestModel).}

\details{
This function is called by the function l() when the argument to l is an ARMA
model. Using l() is usually preferable to calling l.ARMA directly.
l.ARMA calls a compiled program unless compiled=FALSE. The compiled version is much 
faster. 

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
%\keyword{DSE}
\keyword{ts}


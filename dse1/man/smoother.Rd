\name{smoother}
\alias{smoother}
\title{Evaluate a smoother with a TSmodel}
\description{Evaluate a state space model.}
\usage{
    smoother(model, data, compiled=.DSECOMPILED)
}
\arguments{
    \item{model}{
    An object of class `TSestModel' or 'TSmodel' with a model of 
class `non.innov' 'SS' 'TSmodel'. If filter informatin is not provided 
(i.e. in a TSestModel)
then smoother runs the Kalman filter (l.SS) first.}
    \item{data}{A TSdata object.}.
    \item{compiled}{
    If TRUE the compiled version of the code is used. Otherwise the S version is used.}
}
\value{
An object of class TSestModel with an additional element \code{smooth}. 
\code{smooth} is a list of \code{state}, the smoothed state, and 
\code{track}, the smoothed tracking 
error. The result will also contain the element \code{filter} with \code{state}
and \code{track} (which may or may not have been in the original arguement). 
}
\details{
Smoother takes the filtered state \code{estimates$state} and produces a smoothed 
estimate of the state (sometimes called a two sided filter).
}
\seealso{
\code{\link{l}},
\code{\link{l.SS}}
\code{\link{TSmodel}}
\code{\link{TSestModel.object}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
#smoother requires an non-innovations form model
model <- TSmodel(to.SSChol(est.VARX.ls(eg1.DSE.data.diff))) 
smoothed.model <- smoother(model, eg1.DSE.data.diff, compiled=FALSE)
}
%\keyword{DSE}
\keyword{ts}


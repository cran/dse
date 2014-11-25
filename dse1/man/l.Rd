\name{l}
\alias{l}
\alias{l.TSdata}
\alias{l.TSestModel}

\title{Evaluate a TSmodel}
\description{Evaluate a model with data.}
\usage{
    l(obj1, obj2, ...)
    \method{l}{TSdata}(obj1, obj2, ...)
    \method{l}{TSestModel}(obj1, obj2, ...)
}
\arguments{
    \item{obj1}{a TSmodel, TSdata, or TSestModel object.}
    \item{obj2}{a TSmodel or TSdata object.}
    \item{...}{arguments to be passed to other methods.}
}
\value{Usually a TSestModel object is returned. Some methods allow an
  argument \code{result} which specifies that a certain part of the object is
  returned. (For example, the likelihood can be returned. This is useful for
  optimization routines.)}
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
%\keyword{DSE}
\keyword{ts}


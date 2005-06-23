\name{smoother}
\alias{smoother}
\alias{smoother.default}
\alias{smoother.TSestModel}
\title{Evaluate a smoother with a TSmodel}
\description{Evaluate a state space model.}
\usage{
    %%smoother(model, data, compiled=.DSEflags()$COMPILED)
    %%\method{smoother}{default}(model, data, compiled=.DSEflags()$COMPILED)
    %%\method{smoother}{TSestModel}(model, data=TSdata(model),
    %%  compiled=.DSEflags()$COMPILED)
    smoother(model, data, compiled=FALSE)
    \method{smoother}{default}(model, data, compiled=FALSE)
    \method{smoother}{TSestModel}(model, data=TSdata(model), compiled=FALSE)
}
\arguments{
    \item{model}{
    An object of class `TSestModel' or 'TSmodel' with a model of 
class `nonInnov' 'SS' 'TSmodel'. If filter informatin is not provided 
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

Fixed interval smoother for a model as returned by l.SS.
ref. appendix of Shumway and Stoffer,1982, J.of Time Series, 253-264,
       Jazwinski 1970, or Anderson and Moore.
Note: this does not allow the same option as l.SS for calculating over a
   sub-sample. Smoothing is done over the length of the available filter
   data (which will be calculated to the length of the data if not
   supplied). For models with an input smoothing will only be done to the
   length of input data if that is smaller than the available filter data. 
See \link{SS} for details of the model:

       z(t) = Fz(t-1) + Gu(t) + Qe(t)
       y(t) = Hz(t)  + Rw(t)

}
\seealso{
\code{\link{state}},
\code{\link{l}},
\code{\link{SS}}
\code{\link{l.SS}}
\code{\link{TSmodel}}
\code{\link{TSestModel.object}}
}
\examples{
data("eg1.DSE.data.diff", package="dse1")
#smoother requires an non-innovations form model
model <- TSmodel(toSSChol(estVARXls(eg1.DSE.data.diff))) 
smoothed.model <- smoother(model, eg1.DSE.data.diff, compiled=FALSE)
tfplot(state(smoothed.model))
tfplot(state(smoothed.model, filter=TRUE))
#compare
tfplot(state(smoothed.model, smoother=TRUE), state(smoothed.model, filter=TRUE))
}
\concept{DSE}
\keyword{ts}


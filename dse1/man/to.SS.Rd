\name{to.SS}
\alias{to.SS}
\alias{to.SS.TSmodel}
\alias{to.SS.ARMA}
\alias{to.SS.SS}
\alias{to.SS.TSestModel}
\alias{to.SS.augment}
\alias{to.SS.augment.ARMA}
\alias{to.SS.augment.TSestModel}
\alias{to.SS.nested}
\alias{to.SS.nested.SS}
\alias{to.SS.nested.ARMA}
\alias{to.SS.nested.TSestModel}

\title{Convert to State Space Model}
\description{
Convert a model to state space form.
}
\usage{
    to.SS(model, ...)
    \method{to.SS}{TSmodel}(model, n=NULL)
    \method{to.SS}{ARMA}(model, ...)
    \method{to.SS}{SS}(model)
    \method{to.SS}{TSestModel}(model, ...)
    
    to.SS.augment(model, ...)
    \method{to.SS.augment}{ARMA}(model, fuzz=1e-14)
    \method{to.SS.augment}{TSestModel}(model, ...)
    
    to.SS.nested(model, ...)
    \method{to.SS.nested}{ARMA}(model, n=NULL, Aoki=F)
    \method{to.SS.nested}{SS}(model, n=NULL, Aoki=F)
    \method{to.SS.nested}{TSestModel}(model, ...)
    
    }
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{n}{If n is specified then it is used as the state dimension when the 
markov parameter conversion technique is required.}
}
\value{

A state space model in an object of class 'SS' 'TSmodel'.
}
\details{

If the order of the AR polynomial equals or exceeds the MA 
polynomial (and the input polynomial) then the model is converted 
by state augmentation. Otherwise, it is converted by approximating 
the markov coefficients a la Mittnik.  (This may not always work
very well. Compare the results to check.)
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
model <- to.SS(model)
}
%\keyword{DSE}
\keyword{ts}


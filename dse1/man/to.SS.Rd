\name{to.SS}
\alias{to.SS}
\alias{to.SS.ARMA}
\alias{to.SS.SS}
\alias{to.SS.TSestModel}
\alias{to.SSaugment}
\alias{to.SSaugment.ARMA}
\alias{to.SSaugment.TSestModel}
\alias{to.SSnested}
\alias{to.SSnested.SS}
\alias{to.SSnested.ARMA}
\alias{to.SSnested.TSestModel}

\title{Convert to State Space Model}
\description{
Convert a model to state space form.
}
\usage{
    to.SS(model, ...)
    \method{to.SS}{ARMA}(model, ...)
    \method{to.SS}{SS}(model, ...)
    \method{to.SS}{TSestModel}(model, ...)
    
    to.SSaugment(model, ...)
    \method{to.SSaugment}{ARMA}(model, fuzz=1e-14, ...)
    \method{to.SSaugment}{TSestModel}(model, ...)
    
    to.SSnested(model, ...)
    \method{to.SSnested}{ARMA}(model, n=NULL, Aoki=FALSE, ...)
    \method{to.SSnested}{SS}(model, n=NULL, Aoki=FALSE, ...)
    \method{to.SSnested}{TSestModel}(model, ...)
    
    }
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{n}{If n is specified then it is used as the state dimension when the 
      markov parameter conversion technique is required.}
    \item{Aoki}{logical indicating if Aoki's method (which does not work in
       general) should be tried.}
    \item{fuzz}{if the zero lag term of polynomials A and B are within fuzz of
       the identitity matrix then they are not inverted. (i.e. they are assumed
       to be identity.)}
    \item{...}{arguments to be passed to other methods.}
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


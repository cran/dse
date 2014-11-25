\name{to.ARMA}
\alias{to.ARMA}
\alias{to.ARMA.ARMA}
\alias{to.ARMA.SS}
\alias{to.ARMA.TSestModel}

\title{Convert to an ARMA Model}
\description{
Convert a state space model to an ARMA representation. The state is 
eliminated by a method which uses an equivalence that can be demonstrated
by the Cayley Hamilton theorem.
}
\usage{
    to.ARMA(model, ...)
    \method{to.ARMA}{ARMA}(model)
    \method{to.ARMA}{SS}(model, fuzz=1e-10)
    \method{to.ARMA}{TSestModel}(model, ...)
}
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{fuzz}{Parameters closer than fuzz to one or zero are set to 1.0
     or 0.0 respectively}
}
\value{
An object of class 'ARMA' 'TSmodel' containing an ARMA model.
}
\references{See for example Aoki(1989)}
\seealso{
\code{\link{to.SS}}
\code{\link{fix.constants}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
model <- to.ARMA(model)
}
%\keyword{DSE}
\keyword{ts}


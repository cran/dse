\name{toARMA}
\alias{toARMA}
\alias{toARMA.ARMA}
\alias{toARMA.SS}
\alias{toARMA.TSestModel}

\title{Convert to an ARMA Model}
\description{
Convert a state space model to an ARMA representation. The state is 
eliminated by a method which uses an equivalence that can be demonstrated
by the Cayley Hamilton theorem.
}
\usage{
    toARMA(model, ...)
    \method{toARMA}{ARMA}(model, ...)
    \method{toARMA}{SS}(model, fuzz=1e-10, ...)
    \method{toARMA}{TSestModel}(model, ...)
}
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{fuzz}{Parameters closer than fuzz to one or zero are set to 1.0
     or 0.0 respectively}
    \item{...}{arguments to be passed to other methods.}
}
\value{
An object of class 'ARMA' 'TSmodel' containing an ARMA model.
}
\references{See for example Aoki(1989)}
\seealso{
\code{\link{toSS}}
\code{\link{fixConstants}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- toSS(estVARXls(eg1.DSE.data.diff))
model <- toARMA(model)
}
%\keyword{DSE}
\keyword{ts}


\name{checkBalanceMittnik}
\alias{checkBalanceMittnik}
\alias{checkBalanceMittnik.ARMA}
\alias{checkBalanceMittnik.SS}
\alias{checkBalanceMittnik.TSestModel}

\title{Check Balance of a TSmodel}
\description{
Calculate the difference between observability and reachability gramians of the model transformed to Mittnik's form.
}
\usage{
    checkBalanceMittnik(model)
    \method{checkBalanceMittnik}{ARMA}(model)
    \method{checkBalanceMittnik}{SS}(model)
    \method{checkBalanceMittnik}{TSestModel}(model)
}
\arguments{
    \item{model}{An object of class TSmodel.}
}
\value{No value is returned.}
\sideffects{
Differences  between the observability and reachability gramians are printed.}
\details{
Balanced models should have equal observability and reachability gramians.
}
\seealso{
\code{\link{checkBalance}}
\code{\link[dse2]{MittnikReduction}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- toSS(estVARXls(eg1.DSE.data.diff))
    checkBalanceMittnik(model)
}
%\keyword{DSE}
\keyword{ts}


\name{check.balanceMittnik}
\alias{check.balanceMittnik}
\alias{check.balanceMittnik.ARMA}
\alias{check.balanceMittnik.SS}
\alias{check.balanceMittnik.TSestModel}

\title{Check Balance of a TSmodel}
\description{
Calculate the difference between observability and reachability gramians of the model transformed to Mittnik's form.
}
\usage{
    check.balanceMittnik(model)
    \method{check.balanceMittnik}{ARMA}(model)
    \method{check.balanceMittnik}{SS}(model)
    \method{check.balanceMittnik}{TSestModel}(model)
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
\code{\link{check.balance}}
\code{\link[dse2]{reduction.Mittnik}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
    check.balanceMittnik(model)
}
%\keyword{DSE}
\keyword{ts}


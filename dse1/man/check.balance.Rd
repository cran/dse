\name{check.balance}
\alias{check.balance}
\alias{check.balance.SS}
\alias{check.balance.ARMA}
\alias{check.balance.TSestModel}

\title{Check Balance of a TSmodel}
\description{
Calculate the difference between observability and reachability gramians.}
\usage{
    check.balance(model)
    \method{check.balance}{SS}(model)
    \method{check.balance}{ARMA}(model)
    \method{check.balance}{TSestModel}(model)
   }
\arguments{
    \item{model}{A TSmodel object.}
}
\value{No value is returned.}
\sideffects{
Differences  between the observability and reachability gramians are printed.
}
\details{
Balanced models should have equal observability and reachability gramians.
}
\seealso{
    \code{\link{check.balanceMittnik}}
    \code{\link[dse2]{reduction.Mittnik}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- to.SS(est.VARX.ls(eg1.DSE.data.diff))
    check.balance(model)
}
%\keyword{DSE}
\keyword{ts}


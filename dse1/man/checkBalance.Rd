\name{checkBalance}
\alias{checkBalance}
\alias{checkBalance.SS}
\alias{checkBalance.ARMA}
\alias{checkBalance.TSestModel}

\title{Check Balance of a TSmodel}
\description{
Calculate the difference between observability and reachability gramians.}
\usage{
    checkBalance(model)
    \method{checkBalance}{SS}(model)
    \method{checkBalance}{ARMA}(model)
    \method{checkBalance}{TSestModel}(model)
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
    \code{\link{checkBalanceMittnik}}
    \code{\link[dse2]{MittnikReduction}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- toSS(estVARXls(eg1.DSE.data.diff))
    checkBalance(model)
}
%\keyword{DSE}
\keyword{ts}


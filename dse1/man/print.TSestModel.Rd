\name{print.TSestModel}
\alias{print.TSestModel}
\alias{print.SS}
\alias{print.ARMA}

\title{Display TSmodel Arrays}
\description{
Display TSmodel arrays.
}
\usage{
    \method{print}{SS}(x, digits=options()$digits)
    \method{print}{ARMA}(x, digits=options()$digits, L=T, fuzz=1e-10)
    \method{print}{TSestModel}(x, ...)
    }

\arguments{
    \item{x}{An object of class TSmodel or TSestModel.}
    \item{digits}{the number of significant digits}
    \item{L}{
      if T then ARMA model arrays are displayed as a polynomial matrix 
      with L indicating lags. Otherwise, each lag of the array is 
      displayed as a matrix.}
    \item{fuzz}{
      ARMA model polynomial elements with absolute value less than fuzz
      are not displayed (i.e.-as if they are zero)}
    \item{...}{arguments to be passed to other methods.}
}
\value{The object is returned invisibly.}
\sideffects{The model arrays are displayed.}
\note{BUG: digits cannot be controlled for some numbers (e.g.- 1.0 is printed
as 0.9999999999)
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
print(model)
}


\seealso{
    \code{\link{print}}
    \code{\link{summary}}
}
%\keyword{DSE}
\keyword{ts}



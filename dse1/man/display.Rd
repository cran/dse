\name{display}
\alias{display}
\title{Display TSmodel arrays}
\description{
Display TSmodel arrays. (also accessed with the generic function print)
}
\usage{
	display(model, digits=4, L=T, fuzz=1e-10)}
\arguments{

\item{ model}{An object of class TSmodel or TSestModel.}

\item{ digits }{the number of significant digits}
\item{ L}
{
if T then ARMA model arrays are displayed as a polynomial matrix 
with L indicating lags. Otherwise, each lag of the array is 
displayed as a matrix.}
\item{ fuzz}
{
ARMA model polynomial elements with absolute value less than fuzz
are not displayed (i.e.-as if they are zero)}

}
\value{No value is returned}
\sideffects{
The model arrays are displayed.
}
\section{ BUGS} 

digits cannot be controlled for some numbers (e.g.- 1.0 is printed
as 0.9999999999)
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
display(model)
print(model)
}
\keyword{DSE}
\keyword{ts}


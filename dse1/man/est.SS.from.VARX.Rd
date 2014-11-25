\name{est.SS.from.VARX}
\alias{est.SS.from.VARX}

\title{Estimate a state space TSmodel using VAR estimation}
\description{
	Estimate a VAR TSmodel with (optionally) an exogenous input 
	and convert to state space.
}
\usage{
	est.SS.from.VARX(data, ...)}
\arguments{
	\item{data}{
	An object with the structure of an object of class TSdata (see TSdata).}
	\item{...}{See arguements to est.VARX.ls}}

\value{A state space model in an object of class TSestModel.}
\details{This function uses the functions est.VARX.ls and to.SS.}
\section{ REFERENCES} {
	P.D.Gilbert (1993) State Space and ARMA models: An overview of
	estimation and reduction.
}
\seealso{
	\code{\link{to.SS}}
	\code{\link[dse2]{est.SS.Mittnik}}
	\code{\link[dse2]{bft}}
	\code{\link{est.VARX.ls}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <-est.SS.from.VARX(eg1.DSE.data.diff)
}
\keyword{DSE}
\keyword{ts}


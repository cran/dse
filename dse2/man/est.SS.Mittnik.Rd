\name{est.SS.Mittnik}
\alias{est.SS.Mittnik}
\title{Estimate a State Space Model}
\description{Estimate a state space model using Mittnik's markov parameter estimation.
}
\usage{
	est.SS.Mittnik(data, max.lag=6, n=NULL, subtract.means=F, normalize=F)
}
\arguments{
\item{ data}{A TSdata object.}
\item{ max.lag}{The number of markov parameters to estimate.}
\item{ subtract.means}{If T subtract the means from the data before estimation.}
\item{ normalize}{If T normalize the data before estimation.}
}
\value{A state space model in an object of class TSestModel.}
\details{See References}
\section{ REFERENCES} {
	P.D.Gilbert (1993) State Space and ARMA models: An overview of
	estimation and reduction. (Other reference cited therein.)
}
\seealso{
	\code{\link{est.VARX.ls}}
	\code{\link{bft}}
}
\examples{
	if(is.R()) data("egJofF.1dec93.data", package="dse1")
	model <- est.SS.Mittnik(egJofF.1dec93.data)
}
\keyword{DSE}
\keyword{ts}


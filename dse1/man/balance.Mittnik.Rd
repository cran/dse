\name{balance.Mittnik}
\alias{balance.Mittnik}
\alias{balance.Mittnik.svd}

\title{Balance a state space model}
\description{Balance a state space model a la Mittnik.}

\usage{
	balance.Mittnik(model, n=NULL)}
\arguments{
	\item{model}{An TSmodel object.}
	\item{n}
{
n is intended primarily for producing a state space model from the markov
parameters of an ARMA model, but if it is supplied with an SS model the
result will be a model with state dimension n based on the n largest
singular values of the svd of a Hankel matrix of markov parameters generated 
by the original model. If n is not supplied then the singular values are
printed and the program prompts for n.}
}
\value{A state space model balance a la Mittnik in a TSestModel object.}
\details{See references.}
\section{REFERENCES} {
	P.D.Gilbert (1993) State Space and ARMA models: An overview of
	estimation and reduction. (Other reference cited therein)
}
\seealso{
	\code{\link{est.VARX.ls}},
	\code{\link{est.VARX.ar}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- to.SS(TSmodel(est.VARX.ls(eg1.DSE.data.diff)))
	newmodel <-balance.Mittnik(model)
}
\keyword{DSE}
\keyword{ts}


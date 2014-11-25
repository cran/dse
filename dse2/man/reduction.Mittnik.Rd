\name{reduction.Mittnik}
\alias{reduction.Mittnik}
\alias{reduction.Mittnik.from.Hankel}

\title{Balance and Reduce a Model}
\description{
Balance and reduce the state dimension of a state space model a la Mittnik.
}
\usage{
	reduction.Mittnik(model, data=NULL, criterion=NULL, verbose=T,warn=T)
}
\arguments{
\item{model}{An object of class TSmodel or TSestModel.}
\item{data}
{If the supplied model is of class TSestModel and data is not supplied
then it is taken from the model. If the model is of class TSmodel
then data must be supplied.}
}
\value{
A state space model balance a la Mittnik in an object of class TSestModel.
}
\details{
If a state space model is supplied the max. state dimension for the result is
taken from the model. If an ARMA model is supplied then singular values 
will be printed and the program prompts for the max. state dimension.
criterion should be the name of one of the values returned by information.tests,
that is, one of ("port","like","aic","bic","gvc","rice","fpe","taic", 
"tbic","tgvc","trice","tfpe"). If criteria is not specified then the 
program prompts for the state dimension (n) to use for the returned model.
The program requires data to calculate selection criteria. (The program
balance.Mittnik calculates svd criteria only and can be used for reduction
without data.)  
}
\section{ REFERENCES} {
	P.D.Gilbert (1993) State Space and ARMA models: An overview of
	estimation and reduction. (Other reference cited therein)

	S.Mittnik, Multivariate Time Series Analysis With State Space Models, 
	Computers Math Appl. Vol 17, No 8/9, pp1189-1201, 1989.
   
	S.Mittnik, Macroeconomic Forecasting Experience With Balance State Space 
	Models,  International Journal Of Forecasting, Vol 6, pp337-348, 1990.
  
	S.Mittnik, Forecasting With Balanced State Space Representations of 
	Multivariate Distributed Lag Models. J. of Forecasting, Vol.9, 207-218, 1990.  
}
\seealso{
	\code{\link{est.VARX.ls}}
	\code{\link{bft}}
	\code{\link{information.tests}}
}
\examples{
	if(is.R()) data("egJofF.1dec93.data", package="dse1")
	model <- to.SS(est.VARX.ls(egJofF.1dec93.data))
	newmodel <-reduction.Mittnik(model, criterion="taic")
}
\keyword{DSE}
\keyword{ts}


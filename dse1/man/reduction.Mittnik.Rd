\name{reduction.Mittnik}
\alias{reduction.Mittnik}
\alias{reduction.Mittnik.from.Hankel}

\title{Balance and Reduce a Model}
\description{
Balance and reduce the state dimension of a state space model a la Mittnik.
}
\usage{
    reduction.Mittnik(model, data=NULL, criterion=NULL, verbose=TRUE,warn=TRUE)
    reduction.Mittnik.from.Hankel(M, data=NULL, nMax=NULL, 
       criterion=NULL, verbose=TRUE, warn=TRUE, 
       Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
}
\arguments{
    \item{model}{An object of class TSmodel or TSestModel.}
    \item{data}{
      If the supplied model is of class TSestModel and data is not supplied
      then it is taken from the model. If the model is of class TSmodel
      then data must be supplied.}
    \item{criterion}{Criterion to be used for model 
       selection. see \code{information.tests.calculations}.}
    \item{verbose}{logical indicating if information should be printed 
       during estimation.}
    \item{warn}{logical indicating if some warning messages should be suppressed.}
    \item{M}{a matrix. See details.}
    \item{nMax}{integer indicating the state dimension of the largest model
       considered.}
    \item{Spawn}{logical indicating if Splus For loops should be used.}
}
\value{
A state space model balance a la Mittnik in an object of class TSestModel.
}
\details{
\code{reduction.Mittnik} gives nested-balanced state space model using
reduction by svd of the Hankel matrix generated from a model.
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

The function \code{reduction.Mittnik.from.Hankel} is called by
\code{reduction.Mittnik} and typically not by the user, but there are situations
when the former might be called directly. It selects a reduced state 
space model by svd a la Mittnik. Models and several criteria for all 
state dimensions up to the max. state dim. specified are calculated. 
(If nMax is not supplied then svd criteria are printed and the program 
prompts for nMax). The output dimension p is taken from nrow(M).
M is a matrix with p x (m+p)  blocks giving the markov parameters,
that is, the first row of the Hankel matrix. It can be generated from the
model as in the function markov.parms, or from the data, as in the function
est.SS.Mittnik.

\code{data} is necessary only if criteria (AIC,etc) are to be calculated.
}
\references{
    Gilbert, P. D. (1993) State space and ARMA models: An overview of
    the equivalence. Working paper 93-4, Bank of Canada. Available at
    <www.bank-banque-canada.ca/pgilbert>

    Gilbert, P. D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.

    Mittnik, S. (1989), Multivariate Time Series Analysis With State 
    Space Models, Computers Math Appl. Vol 17, No 8/9, pp1189-1201.
   
    Mittnik, S. (1990), Macroeconomic Forecasting Experience With 
    Balance State Space Models,  
    International Journal Of Forecasting, Vol 6, pp337-348.
  
    Mittnik, S. (1990), Forecasting With Balanced State Space 
    Representations of Multivariate Distributed Lag Models. 
    J. of Forecasting, Vol.9, 207-218.  
}
\seealso{
    \code{\link{est.VARX.ls}}
    \code{\link{bft}}
    \code{\link{balance.Mittnik}}
    \code{\link{information.tests}}
    \code{\link[dse1]{information.tests.calculations}}
}
\examples{
    if(is.R()) data("egJofF.1dec93.data", package="dse1")
    model <- to.SS(est.VARX.ls(egJofF.1dec93.data))
    newmodel <-reduction.Mittnik(model, criterion="taic")
}
%\keyword{DSE}
\keyword{ts}


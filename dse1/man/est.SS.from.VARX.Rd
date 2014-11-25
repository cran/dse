\name{est.SS.from.VARX}
\alias{est.SS.from.VARX}

\title{Estimate a state space TSmodel using VAR estimation}
\description{
    Estimate a VAR TSmodel with (optionally) an exogenous input 
    and convert to state space.
}
\usage{
    est.SS.from.VARX(data, warn=T, ...)}
\arguments{
    \item{data}{
    An object with the structure of an object of class TSdata (see TSdata).}
    \item{warn}{Logical indicating if warnings should be printed (T) or
    suppressed (F).}
    \item{...}{See arguements to est.VARX.ls}}

\value{A state space model in an object of class TSestModel.}
\details{This function uses the functions est.VARX.ls and to.SS.}
\references{
    Gilbert, P. D. (1993) State space and ARMA models: An overview of
    the equivalence. Working paper 93-4, Bank of Canada. Available at
    <www.bank-banque-canada.ca/pgilbert>

    Gilbert, P. D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.
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
%\keyword{DSE}
\keyword{ts}


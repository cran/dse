\name{balance.Mittnik}
\alias{balance.Mittnik}
\alias{balance.Mittnik.svd}

\title{Balance a state space model}
\description{Balance a state space model a la Mittnik.}

\usage{
    balance.Mittnik(model, n=NULL)
    balance.Mittnik.svd(M, m, n=NULL)
}
\arguments{
    \item{model}{An TSmodel object.}
    \item{M}{a matrix. See details in \code{reduction.Mittnik}.}
    \item{m}{an integer indicating the number of input series in the model.}
    \item{n}{see details}
}
\value{A state space model in a TSestModel object.}
\details{
\code{balance.Mittnik} calculate a state space model balance a la Mittnik.
n is intended primarily for producing a state space model from the markov
parameters of an ARMA model, but if it is supplied with an SS model the
result will be a model with state dimension n based on the n largest
singular values of the svd of a Hankel matrix of markov parameters generated 
by the original model. If n is not supplied then the singular values are
printed and the program prompts for n. \code{balance.Mittnik} calls
\code{balance.Mittnik.svd}

\code{balance.Mittnik.svd} calculates a nested-balanced state space 
model by svd a la Mittnik. If state dim n is supplied then svd 
criteria are not calculated and the given n is used. Otherwise, 
the singular values are printed and the program prompts for n.
M is a matrix with p x (m+p)  blocks giving the markov parameters,
that is, the first row of the Hankel matrix. It can be generated from the
model as in the function markov.parms, or from the data, as in the function
est.SS.Mittnik.
m is the dimension of input series, which is needed to decompose M.
The output dimension p is taken from nrow(M).

See also \code{reduction.Mittnik} and references.
}
\references{
    Gilbert, P. D. (1993) State space and ARMA models: An overview of
    the equivalence. Working paper 93-4, Bank of Canada. Available at
    <www.bank-banque-canada.ca/pgilbert>

    Gilbert, P.D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.

    S.Mittnik (1989), Multivariate Time Series Analysis With State 
    Space Models, Computers Math Appl. Vol 17, No 8/9, pp1189-1201.
   
    S.Mittnik (1990), Macroeconomic Forecasting Experience With 
    Balance State Space Models,  
    International Journal Of Forecasting, Vol 6, pp337-348.
  
    S.Mittnik (1990), Forecasting With Balanced State Space 
    Representations of Multivariate Distributed Lag Models. 
    J. of Forecasting, Vol.9, 207-218.  
}
\seealso{
    \code{\link{est.VARX.ls}},
    \code{\link{est.VARX.ar}}
    \code{\link{reduction.Mittnik}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- to.SS(TSmodel(est.VARX.ls(eg1.DSE.data.diff)))
    \dontrun{newmodel <-balance.Mittnik(model)}
    # this prints information about singular values and prompts with
    #Enter the number of singular values to use for balanced model:
    # 18 might be a good choice in this example. 
    newmodel <-balance.Mittnik(model, n=18)

}
%\keyword{DSE}
\keyword{ts}



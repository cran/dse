\name{sum.sqerror}
\alias{sum.sqerror}
\title{Calculate sum of squared prediction errors}
\description{

Calculate a weighted sum squared prediction errors for a parameterization.
}
\usage{
    sum.sqerror(parms, model=NULL, data=NULL, error.weights=c(1,1,1,1))}
\arguments{
    \item{para}{A vector of parameters.}
    \item{model}{
    An object of class TSmodel which gives the structure of the model to which
para is applied. \code{model$parms} should be the same length as para.}
    \item{data}{
An object of class TSdata which gives the data with which the model is
to be evaluated.}
    \item{error.weights}{
A vector of weights to be applied to the squared prediction errors.}
}
\value{
The value of the sum squared errors for a prediction horizon given by the 
length of error.weights. Each period ahead is weighted by the corresponding 
weight in error.weights.
}
\details{
This function is primarily for use in parameter optimization,
which requires that an objective function be specified by a vector
of parameters.
}
\seealso{
\code{\link{l}}
\code{\link{l.SS}}
\code{\link{l.ARMA}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
sum.sqerror(1e-10+parms(model), model=TSmodel(model), 
        data=TSdata(model), error.weights=c(1,1,10))
}
%\keyword{DSE}
\keyword{ts}


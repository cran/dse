\name{est.max.like}
\alias{est.max.like}
\alias{est.max.like.TSdata}
\alias{est.max.like.TSmodel}
\alias{est.max.like.TSestModel}
\title{Maximum Likelihood Estimation}
\description{Maximum likelihood estimation.}
\usage{
    est.max.like(obj, ...) 
    \method{est.max.like}{TSmodel}(obj, data, algorithm="optim",
	algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf, hessian=TRUE))
    \method{est.max.like}{TSestModel}(obj, data=TSdata(obj), ...)
    \method{est.max.like}{TSdata}(obj, model, ...) 
}
\arguments{
    \item{obj}{an object of class TSmodel, TSdata or TSestModel}
    \item{data}{a TSdata object.}
    \item{model}{a TSmodel object.}
    \item{algorithm}{the algorithm ('optim', 'nlm' or 'nlmin') to use for maximization.}
    \item{algorithm.args}{arguments for the optimization algorithm.}
    \item{...}{arguments passed on to other methods.}
}
\value{
The value returned is an object of class TSestModel with  additional
elements \code{est$converged}, which is T or F indicating convergence, 
and \code{est$results}. The hessian and gradient in results could potentially 
be used for restarting in the case of non-convergence, but that has not 
yet been implemented.
}
\details{
  Shape is used to specify both the initial parameter values and the model
  structure (the placement of the parameters
  in the various arrays of the TSmodel). Estimation attempts to minimimize the
  negative log likelihood (as returned by \code{l} ) of the given model structure by adjusting the
  parameter values. A TSmodel can also have constant values in some array
  elements, and these are not changed.
}
\warning{
  Maximum likelihood estimation of multivariate time series models tends to be
  problematic, even when a good structure and good starting parameter values
  are known. This is especially true for state space models. Also, it seems that
  in-sample fit is often obtained at the expense of out-of-sample forecasting
  ability. If a prior model structure is not important then the \code{bft}
  estimation method may be preferable.
}
\seealso{
\code{\link{optim}}
\code{\link{nlm}}
\code{\link{est.VARX.ls}}
\code{\link{bft}}
\code{\link{TSmodel}}
\code{\link{l}}
}
\examples{
  true.model <- ARMA(A=c(1, 0.5), B=1)
  est.model <-  est.max.like(true.model,  simulate(true.model))
  summary(est.model)
  est.model
  tfplot(est.model)
}
%\keyword{DSE}
\keyword{ts}


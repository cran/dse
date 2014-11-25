\name{simulate}
\alias{simulate}
\alias{simulate.ARMA}
\alias{simulate.SS}
\alias{simulate.TSestModel}

\title{Simulate a TSmodel}
\description{Simulate a model to produce artificial data.}
\usage{
	simulate(model, input=NULL, input0=NULL,
	    sampleT=100, start. = NULL, freq = NULL, 
	    noise=NULL, noise.model=NULL,sd=1, SIGMA=NULL,
	    rng=NULL, compiled=.DSECOMPILED)
}
\arguments{
  	\item{model}{An object of class TSmodel or TSestModel.}
	\item{input}{Data for the exogenous variable if specified in the model.}
	\item{sampleT}{The length of the sample to simulate.}
	\item{start.}{ start. date for resulting data.}
	\item{freq}{ freq for resulting data.}
	\item{y0, input0}{
Lagged values prior to t=1 for y and u, in reverse order so y0[1,] and
input0[1,]correspond to t=0. These arguments  are not implemented for state 
space models. If not specified initial values are set to zero. }
	\item{noise}{
Noise can be supplied. Otherwise it will be generated. 
If supplied it should be a list as described below under returned value.}
	\item{SIGMA}{
The covariance of the noise process. If this is specified then sd is ignored.
A vector or scalar is treated as a diagonal matrix. For an object of class
TSestModel, if neither SIGMA nor sd are specified, then SIGMA is set to 
the estimated covariance (\code{model$estimates$cov}).}
	\item{sd}{The standard deviation of the noise. This can be a vector.}
	\item{noise.model}{ 
A TSmodel to be used for generating noise (not yet supported by SS methods).}
	\item{rng}{The random number generator and seed can be specified to 
regenerate a simulation.}
	\item{seed}{The seed can be specified to regenerate a simulation. 
(rng is better.)}
	\item{compiled}{
	Specifies the compiled version of the code should be used (instead of 
	the S code version).}
}

\value{
The value returned is an object of class TSdata which can be supplied as an
argument to estimation routines. (See TSdata). In addition to the usual
elements (see the description of a TSdata object)  
there are some additional elements:
   model- the generating model, 
   rng - the initial RNG and seed,
   version - the version of S used (random number generators may vary)
   SIGMA as specified
   sd    as specified
   noise - a list of e, w and w0 - the noise processes. w0 is w for t=0 in
      the state space model and prior lags in ARMA models. For VAR models B has
      no lags so w0 has no effect.
   state - the state variable for state space models.
}
\details{

A state space or ARMA model as described in TSmodel is simulated with pseudo
random noise (The default noise is a normally distributed processes. 
An object of class TSdata is
returned. This can be used as input to estimation algorithms. If start. and 
freq are specified, or if input or \code{noise$w} (in that order)
have time series properties, these are given to the output. 
If \code{noise$w0} is a matrix (rather than a vector) for a state space model 
simulation (as it is for ARMA simulations) then it is set to a
vector of zeros. This provides compatability with VAR models (ARMA
models with no lags in B). In general ARMA and SS simulations will
not produce exactly the same results because it is impossible to
determine necessary transformation of initial conditions and w0.
The rng will be set first if it is specified. 
}
\seealso{
\code{\link{TSmodel}},
\code{\link{TSdata}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <- simulate(model)}
\keyword{DSE}
\keyword{ts}


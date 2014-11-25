\name{MonteCarloSimulations}
\alias{MonteCarloSimulations}
\alias{MonteCarloSimulations.default}
\alias{MonteCarloSimulations.TSmodel}
\alias{MonteCarloSimulations.TSestModel}
\alias{MonteCarloSimulations.EstEval}
\alias{MonteCarloSimulations.monte.carlo.simulation}
\alias{is.monte.carlo.simulation}

\title{Generate simulations}
\description{Run multiple simulations}
\usage{
    is.monte.carlo.simulation(obj)
    MonteCarloSimulations(model, simulation.args=NULL, 
           replications=100, rng=NULL, ...)
    \method{MonteCarloSimulations}{default}(model, simulation.args = NULL, 
 		replications = 100, rng = NULL, quiet = FALSE)
    \method{MonteCarloSimulations}{TSmodel}(model, simulation.args=NULL,
          replications=100, rng=NULL, Spawn=.SPAWN, quiet=F)
    \method{MonteCarloSimulations}{TSestModel}(model, simulation.args=NULL, 
           replications=100, rng=NULL, ...)
    \method{MonteCarloSimulations}{EstEval}(model,...)
    \method{MonteCarloSimulations}{monte.carlo.simulation}(model,...)
}
\arguments{
    \item{model}{A model with a simulate method (e.g. a TSmodel).}
    \item{simulation.args,}{
    A list of arguments in addition to model which are passed to 
    simulate.}
    \item{replications}{The number of simulations.}
    \item{rng}{The RNG and starting seed.}
    \item{Spawn}{If T "For" loops are used in Splus.}
    \item{quiet}{logical indicating if printing and many warning messages should
       be suppressed.}
    \item{obj}{an object.}
    \item{...}{arguments passed to other methods.}
}
\value{A list of simulations.}
\details{
This function runs many simulations using \code{simulate}.
Often it not be necessary to do this since the seed can be used to
reproduce the sample and many functions for testing estimation methods, etc.,
will produce samples as they proceed. This function is useful for verification
and for looking at the stochastic properties of the output of a model.
If \code{model} is an object of class \code{EstEval} or 
\code{simulation}
then the model and the seed!!! are extracted so the same sample will be 
generated. The default method expects the result of \code{simulate(model)} to be
a matrix.
There is a \code{tfplot} method (time series plots of the simulations) and a
\code{distribution} method for the result. The latter plots kernel estimates 
of the distribution of the simulations at specified periods.
}
\seealso{
\code{\link{simulate}}
\code{\link{EstEval}}
\code{\link{distribution}}
\code{\link{forecastCov.wrt.true}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <-  MonteCarloSimulations(model, simulation.args=list(sampleT=100))
tfplot(z)
distribution(z)
}
%\keyword{DSE}
\keyword{ts}


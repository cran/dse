\name{distribution.MonteCarloSimulations}
\alias{distribution.MonteCarloSimulations}
\title{Generate distribution plots of Monte Carlo simulations}
\description{Generate distribution plots of Monte Carlo simulations.}
\usage{
    \method{distribution}{MonteCarloSimulations}(obj,
        series=seq(dim(obj$simulations)[2]),
        x.sections=TRUE, periods=1:3, graphs.per.page=5)
}
\arguments{
    \item{obj}{The result of monte.carlo.simulation.}
    \item{series}{The series which should be plotted. The default gives all series.}
    \item{x.sections}{If TRUE then kernel density estimates are plotted 
       for periods indicated by periods. If FALSE then a time series plots 
       of the mean and estimates 1  and 2 standard deviations from the 
       mean. Periods is ignored if x.sections is FALSE.}
    \item{periods}{
       The periods at which the distribution should be calculated and plotted.
       The default gives the first three.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
}
\value{None}
\details{
Kernel estimates of the densities (series by series, not joint densities) are estimated using ksmooth (if available) or density (if available) to produces density plots.
Output graphics can be paused between pages by setting par(ask=TRUE).
}
\seealso{
\code{\link{tfplot.MonteCarloSimulations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <-  MonteCarloSimulations(model)
distribution(z)
}
%\keyword{DSE}
\keyword{ts}


\name{tfplot.monte.carlo.simulations}
\alias{tfplot.monte.carlo.simulations}
\title{Generate plots of Monte Carlo simulations}
\description{Generate plots of Monte Carlo simulations.}
\usage{
    \method{tfplot}{monte.carlo.simulations}(x, start.=NULL, end.=NULL,
    select.series=seq((dim(x$simulations)[2])), 
    select.simulations=seq(dim(x$simulations)[3]),
    graphs.per.page=5, mar=par()$mar)
}
\arguments{
    \item{obj}{The result of monte.carlo.simulation.}
    \item{start.}{
    The starting period for plots. (The default NULL gives the start of 
    the simulations.)}
    \item{end.}{
    The ending period for plots.(The default NULL gives the end of 
    the simulations.)}
    \item{select.series}{
    The series which should be plotted. The default NULL gives all series.}
    \item{select.simulations}{
    Vector of integers indicating the simulations which should be plotted. 
    The default plots all simulations.}
    \item{graphs.per.page}{The number of graphs to put on a page.}
    \item{mar}{Plot margins (see \code{par}).}
}
\value{None}
\details{
This function produces plots of the simulated series.
Output graphics can be paused between pages by setting par(ask=T).
}
\seealso{
\code{\link{distribution.monte.carlo.simulations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <-  monte.carlo.simulations(model)
tfplot(z)
}
%\keyword{DSE}
\keyword{ts}


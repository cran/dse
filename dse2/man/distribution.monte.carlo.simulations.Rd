\name{distribution.monte.carlo.simulations}
\alias{distribution.monte.carlo.simulations}
\title{Generate distribution plots of Monte Carlo simulations}
\description{Generate distribution plots of Monte Carlo simulations.}
\usage{
	distribution(obj, select.series=seq(dim(obj$simulations)[2]),
	     x.sections=T, periods=1:3, graphs.per.page=5)
}
\arguments{
	\item{obj}{The result of monte.carlo.simulation.}
	\item{select.series}{The series which should be plotted. The default gives all series.}
	\item{x.sections}{If T then kernel density estimates are plotted 
	for periods indicated by periods. If F then a time series plots 
	of the mean and estimates 1  and 2 standard deviations from the 
	mean. Periods is ignored if x.sections is F.}
	\item{periods}{
	The periods at which the distribution should be calculated and plotted.
	The default gives the first three.}
	\item{graphs.per.page}{The number of graphs to put on a page.}
}
\value{None}
\details{
Kernel estimates of the densities (series by series, not joint densities) are estimated using ksmooth (if available) or density (if available) to produces density plots.
Output graphics can be paused between pages by setting par(ask=T).
}
\seealso{
\code{\link{tfplot.monte.carlo.simulations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
z <-  monte.carlo.simulations(model)
distribution(z)
}
\keyword{DSE}
\keyword{ts}


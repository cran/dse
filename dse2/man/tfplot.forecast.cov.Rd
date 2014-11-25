\name{tfplot.forecast.cov}
\alias{tfplot.forecast.cov}
\title{Plots of Forecast Variance}
\description{Generate plots of forecast variance calculated by forecast.cov.}
\usage{
    \method{tfplot}{forecast.cov}(x,
       select.series = 1:dim(x$forecast.cov[[1]])[2], 
       select.cov = 1:length(x$forecast.cov), select.true = T, 
       select.zero = T, select.trend = T, y.limit = NULL, line.labels = F, 
       lty = NULL, Legend = NULL, Title = NULL,
       graphs.per.page = 5, mar=par()$mar,  ...)
}
\arguments{
    \item{obj}{The result of forecast.cov.}
    \item{select.series}{
      The series which should be plotted. (The default NULL gives 
      the end of the simulations.)}
    \item{select.cov}{
      For the case of multiple models select the covariance to be plotted.}
    \item{graphs.per.page}{The number of graphs to put on a page.}
    \item{mar}{Plot margins (see \code{par}).}
}
\value{None}
\details{
    This function produces plots of the variance at different horizons.
Output graphics can be paused between pages by setting par(ask=T).
}
\seealso{
    \code{\link{plot}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    z <- forecast.cov(model, eg1.DSE.data.diff)
    tfplot(z)
}
%\keyword{DSE}
\keyword{ts}


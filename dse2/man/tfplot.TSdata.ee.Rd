\name{tfplot.TSdata.ee}
\alias{tfplot.TSmodel.ee}
\alias{tfplot.TSestModel.ee}
\alias{tfplot.featherForecasts}
\alias{tfplot.forecast}
\alias{tfplot.forecastCov.estimatorsWRTdata}
\alias{tfplot.horizonForecasts}
\alias{tfplot.multi.model.horizonForecasts}

\title{Specific Methods for tfplot}
\description{See the generic function description.}
\usage{
    \method{tfplot}{TSmodel.ee}(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL)
    \method{tfplot}{TSestModel.ee}(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL)
    \method{tfplot}{featherForecasts}(x, start.=NULL, end.=NULL, 
      series=NULL, graphs.per.page=5, reset.screen=T, mar=par()$mar)
    \method{tfplot}{forecast}(x, start.=NULL, end.=NULL,
        series = seq(length=nseriesOutput(x$data)),
        names = seriesNamesOutput(x$data), mar=par()$mar)
    \method{tfplot}{forecastCov.estimatorsWRTdata}(x, 
        series=1:dim(x$forecastCov[[1]])[2], 
        select.cov=1:length(x$forecastCov),
        select.zero=T, select.trend=T, lty=NULL, mar=par()$mar,  ...)
    \method{tfplot}{horizonForecasts}(x, start.=NULL, end.=NULL,
         series=NULL, names=seriesNamesOutput(x$data), mar=par()$mar)
    \method{tfplot}{multi.model.horizonForecasts}(x, start.=NULL, end.=NULL,
         series=NULL)
}
\arguments{
    \item{x}{an object for which a tfplot is to be produced.}
    \item{start.}{see \code{tfplot}.}
    \item{end.}{see \code{tfplot}.}
    \item{names}{vector of strings for y axis labelling (ylab).}
    \item{...}{object or objects to be plotted.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{reset.screen}{logical indicating if the plot window should be cleared
       before starting.}

    \item{series}{integer or string indicating
      the series which should be plotted.}
    \item{select.cov}{logical indicating that
      for the case of multiple models select the covariance to be plotted.}
    \item{select.zero}{logical indicating that results from a forecast of 
      zero should be plotted.}
    \item{select.trend}{logical indicating that results from a forecast of 
      trend should be plotted.}
    \item{mar}{plot margins. See \code{par}.}
    \item{graph.args}{list of graphics arguments eventually passed 
       to \code{plot}. See \code{par}.}
    \item{lty}{optional. see plot.}
    \item{criterion}{criterion which should be used to extract something from
         the object which will then be plotted. See \code{EstEval}.}
    \item{criterion.args}{arguments to be passed to code{criterion}.}
}

\seealso{
    \code{\link{tfplot}}
    \code{\link{EstEval}} }

%\keyword{DSE}
\keyword{ts}


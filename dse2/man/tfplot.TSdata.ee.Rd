\name{tfplot.TSdata.ee}
\alias{tfplot.TSmodelEstEval}
\alias{tfplot.TSestModelEstEval}
\alias{tfplot.featherForecasts}
\alias{tfplot.forecast}
\alias{tfplot.EstEval}
\alias{tfplot.horizonForecasts}
\alias{tfplot.multiModelHorizonForecasts}

\title{Specific Methods for tfplot}
\description{See the generic function description.}
\usage{
    \method{tfplot}{TSmodelEstEval}(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL)
    \method{tfplot}{TSestModelEstEval}(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL)
    \method{tfplot}{featherForecasts}(x, start.=NULL, end.=NULL, 
        series=seq(nseries(x)), 
        Title="Predictions (dotted) and actual data (solid)", 
        ylab=seriesNamesOutput(x),
        graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
    \method{tfplot}{forecast}(x, start.=NULL, end.=NULL,
        series = seq(length=nseriesOutput(x$data)),
	Title="Predictions (dotted) and actual data (solid)",
        ylab = seriesNamesOutput(x$data), 
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
    \method{tfplot}{EstEval}(x, start.=NULL, end.=NULL,
        series = seq(length=nseriesOutput(x$truth)),
	Title="Predictions (dotted) and actual data (solid)",
        ylab = seriesNamesOutput(x$truth), 
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
    \method{tfplot}{horizonForecasts}(x, start.=NULL, end.=NULL,
         series=seq(length=nseriesOutput(x$data)),
	 Title="Predictions (dotted) and actual data (solid)", 
	 ylab=seriesNamesOutput(x$data), 
	 graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
    \method{tfplot}{multiModelHorizonForecasts}(x, start.=NULL, end.=NULL,
         series=NULL)
}
\arguments{
    \item{x}{an object for which a tfplot is to be produced.}
    \item{start.}{see \code{tfplot}.}
    \item{end.}{see \code{tfplot}.}
    \item{Title}{string of characters to use for title.}
    \item{ylab}{vector of strings for y axis labelling.}
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
    \code{\link[tframe]{tfplot}}
    \code{\link{EstEval}} }

%\keyword{DSE}
\keyword{ts}


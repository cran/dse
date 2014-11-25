\name{tfplot.TSdata}
\alias{tfplot.TSmodel.ee}
\alias{tfplot.TSestModel.ee}
\alias{tfplot.feather.forecasts}
\alias{tfplot.forecast}
\alias{tfplot.forecast.cov.estimators.wrt.data}
\alias{tfplot.horizon.forecasts}
\alias{tfplot.multi.model.horizon.forecasts}
\alias{tfplot.parms.ee}
\alias{tfplot.roots.ee}
\alias{plot.roots.ee}

\title{Specific Methods for tfplot}
\description{See the generic function description.}
\usage{
    \method{tfplot}{TSmodel.ee}(x, graph.args=NULL,
                       criterion ="parms", criterion.args=NULL)
    \method{tfplot}{TSestModel.ee}(x, graph.args=NULL,
                       criterion ="parms", criterion.args=NULL)
    \method{tfplot}{feather.forecasts}(x, start.=NULL, end.=NULL, 
      select.series=NULL, graphs.per.page=5, reset.screen=T, mar=par()$mar)
    \method{tfplot}{forecast}(x, start.=NULL, end.=NULL,
        select.series = seq(length=output.dimension(x$data)),
        names = output.series.names(x$data), mar=par()$mar)
    \method{tfplot}{forecast.cov.estimators.wrt.data}(x, 
        select.series=1:dim(x$forecast.cov[[1]])[2], 
        select.cov=1:length(x$forecast.cov),
        select.zero=T, select.trend=T, lty=NULL, mar=par()$mar,  ...)
    \method{tfplot}{horizon.forecasts}(x, start.=NULL, end.=NULL,
         select.series=NULL, names=output.series.names(x$data), mar=par()$mar)
    \method{tfplot}{multi.model.horizon.forecasts}(x, start.=NULL, end.=NULL,
         select.series=NULL)
    \method{tfplot}{parms.ee}(x, cum=T, norm=F, bounds=T, invert=F, Sort=F,
	graphs.per.page = 5)
    \method{plot}{roots.ee}(x, complex.plane=T, cum=T, norm=F, bounds=T,
      transform=NULL, invert=F, Sort=T)
    \method{tfplot}{roots.ee}(x, ...)
}

\seealso{ \code{\link{tfplot}} }

%\keyword{DSE}
\keyword{ts}


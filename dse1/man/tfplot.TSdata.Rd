\name{tfplot.TSdata}
\alias{tfplot.TSdata}
\alias{tfplot.TSmodel}
\alias{tfplot.TSestModel}

\title{Tfplot Specific Methods}
\description{See the generic function description.}
\usage{
    \method{tfplot}{TSdata}(..., start. = NULL, end. = NULL, Title = "", 
        reset.screen=TRUE,
        select.inputs = seq(length= nseriesInput(data)),
        select.outputs= seq(length=nseriesOutput(data)),
	mar=par()$mar, 
        graphs.per.page=5, ylab=NULL)
    \method{tfplot}{TSestModel}(..., start.=NULL,end.=NULL,Title=NULL, 
        reset.screen=T, select.inputs=NULL, select.outputs=NULL,
        graphs.per.page=5, mar=par()$mar)
    }
\arguments{
    \item{...}{object or objects to be plotted.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{reset.screen}{logical indicating if the plot window should be cleared
       before starting.}
    \item{Title}{string to use for plot title (passed to plot).}
    \item{ylab}{string to use for y label (passed to plot).}
    \item{start.}{start of plot. (passed to tfwindow)}
    \item{end.}{end of plot. (passed to tfwindow)}
    \item{select.inputs}{series to be plotted. (passed to select.series)}
    \item{select.outputs}{series to be plotted. (passed to select.series)}
    \item{mar}{margins passed to plot. See par.)}
}

\seealso{ \code{\link{tfplot}} }

%\keyword{DSE}
\keyword{ts}


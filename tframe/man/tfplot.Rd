\name{tfplot}
\alias{tfplot}
\alias{tfplot.default}
\alias{tfOnePlot}

\title{Plot Tframed Objects}
\description{
    Plot tframe or tframed objects. 
}
\usage{
    tfplot(x, ...)

    \method{tfplot}{default}(x, ..., tf=tfspan(x, ...), start=tfstart(tf), end=tfend(tf),
       series=seq(nseries(x)), Title=NULL, xlab=NULL, ylab=seriesNames(x), 
       graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
    tfOnePlot(x, xlab=NULL, ylab=NULL, 
              tf=tframe(x), start=tfstart(tf), end=tfend(tf),...)
}

\arguments{
    \item{x}{a tframe or tframed object to plot.}
    \item{...}{any additional tframed objects for the same plot.}
    \item{start}{start of plot. (passed to tfwindow)}
    \item{end}{end of plot. (passed to tfwindow)}
    \item{tf}{a tframe or tframed object which can be used to specify start and
    end.}
    \item{series}{series to be plotted. (passed to selectSeries)}
    \item{Title}{string to use for plot title.}
    \item{xlab}{string to use for x label (passed to plot).}
    \item{ylab}{string to use for y label (passed to plot).}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{mar}{margins passed to plot. See par.)}
    \item{reset.screen}{logical indicating if the plot window should be cleared
       before starting. If this is not TRUE then mar values will have no effect.}
}
\value{None.}
\details{
    In many cases these are the same as plot methods, however, tfplot puts
    different series in the object \code{x} in different plot panels, whereas
    \code{plot} usually puts them in the same panel. For this reason,
    \code{tfplot} tends to work better when the scale of the different series
    are very different. If additional objects are supplied, then they should
    each have the same number of series as \code{x} and all corresponding series
    will be plotted in the same panel.
    
    \code{tfplot} provides an alternate generic
    mechanism for plotting time series data.
    New classes of time series may define
    there own \code{tfplot} (and \code{plot}) methods.
    
    The start and end arguments to tfplot 
    determine the start and end of the plot. The argument tf is an alternate way
    to specify the start and end. It is ignored if start and end are specified.
}
\sideffects{An object is plotted.}
\seealso{
    \code{\link{tfprint}}
    \code{\link{tframe}}
    \code{\link{tframed}}
    \code{\link{print}}
    \code{\link{plot}}
}
\examples{
    tfplot(ts(rnorm(100), start=c(1982,1), frequency=12))
    tfplot(ts(rnorm(100), start=c(1982,1), frequency=12), start=c(1985,6))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{hplot}



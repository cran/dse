\name{tfplot}
\alias{tfplot}
\alias{tfplot.default}
\alias{tfOnePlot}

\title{Plot Tframed Objects}
\description{
    Plot tframe or tframed objects. 
}
\usage{
    tfplot(obj, ...)

    \method{tfplot}{default}(..., xlab=NULL, ylab=NULL, graphs.per.page = 5,
                         start.=NULL, end.=NULL,
             series=seq(nseries(list(...)[[1]])), mar=par()$mar)

    tfOnePlot(obj, xlab=NULL, ylab=NULL, start.=NULL, end.=NULL, ...)
}

\arguments{
    \item{obj}{a tframe or tframed object to plot.}
    \item{...}{any additional tframe or tframed objects for the same plot.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{xlab}{string to use for x label (passed to plot).}
    \item{ylab}{string to use for y label (passed to plot).}
    \item{start.}{start of plot. (passed to tfwindow)}
    \item{end.}{end of plot. (passed to tfwindow)}
    \item{series}{series to be plotted. (passed to select.series)}
    \item{mar}{margins passed to plot. See par.)}
}
\value{None.}
\details{
    In many cases these are the same as plot methods. 
    However, tfplot tries to provide an alternate generic
    mechanism that is consistent with the tframe view of the data.
    This may not always be the preferred print and plot method. 
    Also, new classes of time series may define
    there own plot methods in ways which use a different logic
    from the tframe library. Thus tfplot provides a
    way to program functions which use methods consistent with the tframe
    library logic.
    
    Note that the start and end arguments to tfplot are "start." and "end.".
    In some situations truncation works and these can be abreviated, but
    often it does not work and the error message is opague.
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
    tfplot(ts(rnorm(100), start=c(1982,1), frequency=12), start.=c(1985,6))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{hplot}



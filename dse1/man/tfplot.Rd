\name{tfplot}
\alias{tfplot}
\title{Time Series Plots}
\description{
Plot actual data and one-step ahead predictions for one or multiple models.
}
\usage{
tfplot(..., start.=start(TSdata(model)),end.=end(TSdata(model)),
    select.inputs=NULL, select.outputs=NULL, graphs.per.page=5)
}
\arguments{
\item{model}{Objects of class TSestModel or TSdata.}
\item{start.}{The start of the plot.}
\item{end.}{The end of the plot.}
\item{select.inputs}
{
Selected inputs to plot. The default plots all. If plotting is to be
suppressed then this argument should be given the value 0.}
\item{ select.outputs}
{
Selected outputs to plot. The default plots all. If plotting is to be
suppressed then this argument should be given the value 0.}
\item{ graphs.per.page}{Sets the number of graphs to plot per page.}
}
\value{No returned value.}
\details{
If the arguments are of class TSdata then the input and output data are
plotted (subject to select.inputs and select.outputs).
If the arguments are of class TSestModel then model one step ahead predictions
are also plotted.
Use dev.ask(T) to pause plotting between pages.
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
tfplot(eg1.DSE.data.diff)
model <- est.VARX.ls(eg1.DSE.data.diff)
tfplot(model)
}
\keyword{DSE}
\keyword{ts}


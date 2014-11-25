\name{phasePlots}
\alias{phasePlots}
\title{Calculate Phase Plots}
\description{Calculate phase plots}
\usage{
    phasePlots(data, max.lag=1,diff=FALSE)
}
\arguments{
    \item{data}{A matrix, time series matrix, or an object of class TSdata.}
    \item{max.lag}{The maximum number of shifts to plot}
    \item{diff}{If TRUE the data is plotted against the difference with lagged values.}
}
\value{None}
\sideffects{
A plot of (the phase space) the data against (differenced) lagged values is produced.
}
\details{
Non-linearities may show up as a non-linear surface, but this is
a projection so, for example, a spherical space would not show up.
Some sort of cross-section window would show this but require even more plots.
A good statistical test would be better!
}
\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
phasePlots(egJofF.1dec93.data)
}
%\keyword{DSE}
\keyword{ts}


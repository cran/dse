\name{shock.decomposition}
\alias{shock.decomposition}
\title{Shock Decomposition}
\description{Graphs of the effect of shocks are plotted.}
\usage{
    shock.decomposition(model, horizon=30, shock=rep(1,horizon))
}
\arguments{
    \item{model}{An object of class TSmodel or TSestModel.}
    \item{horizon}{The number of periods for which to calculate the effect of shocks.}
}
\value{None}
\sideffects{
Graphs of the effect of shocks are plotted. 
}
\details{All output data is set to zero and then each output in turn is switched to 
a value of 1.0 for all periods. 
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
shock.decomposition(model)
}
%\keyword{DSE}
\keyword{ts}


\name{shockDecomposition}
\alias{shockDecomposition}
\title{Shock Decomposition}
\description{Graphs of the effect of shocks are plotted.}
\usage{
    shockDecomposition(model, horizon=30, shock=rep(1,horizon))
}
\arguments{
    \item{model}{An object of class TSmodel or TSestModel.}
    \item{horizon}{The number of periods for which to calculate the effect 
        of shocks.}
    \item{shock}{data to be used model output. See details.}
}
\value{None}
\sideffects{
Graphs of the effect of shocks are plotted. 
}
\details{All output data is set to zero and then each output in turn is switched to 
a value of shock (default 1.0) for all periods. 
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- estVARXls(eg1.DSE.data.diff)
shockDecomposition(model)
}
%\keyword{DSE}
\keyword{ts}


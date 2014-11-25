\name{est.SS.Mittnik}
\alias{est.SS.Mittnik}
\title{Estimate a State Space Model}
\description{Estimate a state space model using Mittnik's markov parameter estimation.
}
\usage{
    est.SS.Mittnik(data, max.lag=6, n=NULL, subtract.means=F, normalize=F)
}
\arguments{
    \item{data}{A TSdata object.}
    \item{max.lag}{The number of markov parameters to estimate.}
    \item{n}{The state dimension.}
    \item{subtract.means}{If T subtract the means from the data before estimation.}
    \item{normalize}{If T normalize the data before estimation.}
}
\value{A state space model in an object of class TSestModel.}
\details{See reduction.Mittnik and references cited there. If the state
dimension is not specified then the singular values of the Hankel matrix are
printed and the user is prompted for the state dimension.
}
\references{
    See references for reduction.Mittnik. 
}
\seealso{
    \code{\link{reduction.Mittnik}}
    \code{\link{est.VARX.ls}}
    \code{\link{bft}}
}
\examples{
    if(is.R()) data("egJofF.1dec93.data", package="dse1")
    \dontrun{model <- est.SS.Mittnik(egJofF.1dec93.data)}
    # this prints information about singular values and prompts with
    #Enter the number of singular values to use for balanced model:
    # the choice is difficult in this example. 
    model <- est.SS.Mittnik(egJofF.1dec93.data, n=3)
}
%\keyword{DSE}
\keyword{ts}


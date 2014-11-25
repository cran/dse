\name{est.black.box}
\alias{est.black.box}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    est.black.box(data,...)
}
\arguments{
    \item{data}{Data in an object of class TSdata.}
    \item{...}{
    Optional arguments depent on the function which is eventually called.}
}
\value{A state space model in an object of class TSestModel.}
\details{
The function makes a call the most promising procedure currently available.
These tend to have names like est.black.box1, est.black.box2, est.black.box4.
This is an active area of ongoing research and so the actual routine called
will probably change with new versions.
}
\examples{
if(is.R()) data("egJofF.1dec93.data", package="dse1")
goodmodel <- est.black.box(egJofF.1dec93.data)
}
%\keyword{DSE}
\keyword{ts}


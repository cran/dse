\name{estBlackBox}
\alias{estBlackBox}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    estBlackBox(data,...)
}
\arguments{
    \item{data}{Data in an object of class TSdata.}
    \item{...}{
    Optional arguments depent on the function which is eventually called.}
}
\value{A state space model in an object of class TSestModel.}
\details{
The function makes a call the most promising procedure currently available.
These tend to have names like estBlackBox1, estBlackBox2, estBlackBox4.
This is an active area of ongoing research and so the actual routine called
will probably change with new versions.
}
\examples{
data("egJofF.1dec93.data", package="dse")
goodmodel <- estBlackBox(egJofF.1dec93.data)
}
\concept{DSE}
\keyword{ts}


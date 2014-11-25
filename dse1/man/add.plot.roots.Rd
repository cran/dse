\name{add.plot.roots}
\alias{add.plot.roots}
\title{Add Model Roots to a plot}
\description{

Calculate and plot roots of a model.
}
\usage{
    add.plot.roots(v, pch='*', fuzz=0)
}
\arguments{

    \item{v}{An object containing a TSmodel.}
    \item{pch}{Character to use for plotting.}
    \item{fuzz}{If non-zero then roots within fuzz distance are considered equal.}
}
\sideffects{
The roots are addeded to an existing plot.
}
\value{
The eigenvalues of the state transition matrix or the inverse of the roots of the
determinant of the AR polynomial are returned invisibly.
}
\seealso{
    \code{\link{plot.roots}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
plot.roots(model)
add.plot.roots(to.SS(model))
}
%\keyword{DSE}
\keyword{ts}


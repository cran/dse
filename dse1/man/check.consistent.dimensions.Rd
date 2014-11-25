\name{check.consistent.dimensions}
\alias{check.consistent.dimensions}
\alias{check.consistent.dimensions.default}
\alias{check.consistent.dimensions.ARMA}
\alias{check.consistent.dimensions.SS}
\alias{check.consistent.dimensions.TSdata}
\alias{check.consistent.dimensions.TSestModel}

\title{Check Consistent Dimensions}
\description{Check that dimensions of a model and data agree.}
\usage{
	check.consistent.dimensions(model,data)}
\arguments {
	\item{model}{A TSmodel}
	\item{data}{TSdata}
}
\value{logical}
\details{
	Check that dimensions of a model and data agree.
}


\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.VARX.ls(eg1.DSE.data.diff)
	check.consistent.dimensions(model)
}
\keyword{DSE}
\keyword{ts}


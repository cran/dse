\name{standardize}
\alias{standardize}
\title{standardize data}
\description{Convert data to mean zero, standard deviation one. This does not
remove cross correlations between series.}
\usage{standardize(ser)}
\arguments{
    \item{ser}{A vector time matrix.}
    }
\value{A time series matrix.}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- standardize(outputData(eg1.DSE.data.diff))
}
\concept{DSE}
\keyword{ts}


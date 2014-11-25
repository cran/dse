\name{diff.log}
\alias{diff.log}
\title{Calculate the difference of log data}
\description{Calculate the difference from lag periods prior for log of data.}
\usage{
	diff.log(series, lag=1, base = 2.71828182845905)}
\arguments{
\item{ series}{A time series.}
\item{ lag}{The difference is calculated relative to lag periods prior.}
}
\value{
A time series of the difference relative to lag periods prior for the log of 
the data. lag data points are lost from the beginning of the series. Negative
values will result in NAs.
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
z <- diff.log(output.data(eg1.DSE.data))
}
\keyword{DSE}
\keyword{ts}


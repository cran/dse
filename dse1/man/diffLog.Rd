\name{diffLog}
\alias{diffLog}
\title{Calculate the difference of log data}
\description{Calculate the difference from lag periods prior for log of data.}
\usage{
    diffLog(x, lag=1, base = exp(1), ...)}
\arguments{
    \item{x}{A time series.}
    \item{lag}{The difference is calculated relative to lag periods prior.}
    \item{base}{Base to use when calculating logrithms.}
    \item{...}{(further arguments, currently disregarded).}
}
\value{
A time series of the difference relative to lag periods prior for the log of 
the data. lag data points are lost from the beginning of the series. Negative
values will result in NAs.
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
z <- diffLog(outputData(eg1.DSE.data))
}
\concept{DSE}
\keyword{ts}


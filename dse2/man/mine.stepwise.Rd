\name{mine.stepwise}
\alias{mine.stepwise}
\title{Mine Stepwise}
\description{mine.stepwise}
\usage{
	mine.stepwise(data, essential.data=1,
      method="efroymson", f.crit=2, intercept=T,
      subtract.means=F,  standardize=F, 
      lags.in=6, lags.out=6, trend=F, plot.=T)
}
\arguments {
	\item{data}{TSdata}
	\item{essential.data}{An integer vector indication important data.}
	\item{method}{method to pass to stepwise.}
	\item{f.crit}{See details.}
	\item{intercept}{See details.}
	\item{subtract.means}{See details.}
	\item{standardize}{See details.}
	\item{lags.in}{See details.}
	\item{lags.out}{See details.}
	\item{trend}{See details.}
	\item{plot.}{See details.}
}
\value{x}
\details{This documentation is out of date.
	Data should be of class TSdata. essential.data must have length 1.
	standardize and subtract means ...
	The result is a list with the results of stepwise,..., and
	several vectors indicating information about the columns of the matrix
	passed to stepwise: 
	io.indicator indicating an input (F) or output (T) variable
	v.indicator  indicating which series
	lag.indicator indicating the lag
	s.input.indicator and s.output.indicator are logic matrices
	length(stepwise$rss) by m and p respectively indicating if
	a series is included for each element of rss.
}

\keyword{DSE}
\keyword{ts}


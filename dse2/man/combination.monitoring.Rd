\name{combination.monitoring}
\alias{combination.monitoring}
\alias{combine.and.forecast}
\alias{reconstruct.combined.forecast}
\alias{construct.data.to.override.horizon}
\alias{get.overriding.data}
\alias{restrict.overriding.data}

\title{Combination Monitoring}
\description{Automatic monitoring with e-mail of results}
\usage{
	combination.monitoring <- function(model, data.names,
	   previous.data=NULL,
	   overriding.data.names=NULL, 
	   restrict.overriding.data=T, overriding.horizon=0,
	   mail.list=NULL,
	   error.mail.list=NULL,
	   message.title="Combination Monitoring",
	   message.subject="Combination Monitoring",
	   message.footnote=NULL,
	   show.start= c(0,-3),
	   show.end  = c(0,12),    
	   report.variables=series.names(data.names),
	   data.sub.heading=NULL,
	   data.tag=" ",
	   future.input.data.tag="p",
	   overriding.data.tag="m",
	   overlapping.period.forecast.tag="g",
	   forecast.tag="f",
	   run.again=F,
	   save.as=NULL)
}
\arguments {
	\item{x}{see details below}
}
\value
	Invisibly return latest data for next comparison. This function is
	run mainly for its side effects.}
\details{
	This function allows for 
	   -combining forecasts (ie. monitoring or other forecast data)
	   -input (policy) projections
	If these feature are not need see simple.monitoring.
	mail.list and error.mail.list should be single strings (not vectors)
	but the string can contain multiple user ids for mail
	If overriding.data.names=NULL then no overriding data is used.
	report.variables indicates output variables which are reported. If NULL,
	then all outputs are reported.
	show.end is min of this and overriding.data$input if needed.

	This functions allow for the use of over-riding data
	which may come from other forecasts or monitoring and 
	can be used to augment (and replace) actual data.
	Also, input (policy) variable forecasts can be used.<
	NB. The combination is not in the sense of averaging
	together forecasts.    

}
\seealso{
	\code{\link{simple.monitoring}}
	\code{\link{combination.monitor.function.tests}} 
}

\keyword{DSE}
\keyword{ts}


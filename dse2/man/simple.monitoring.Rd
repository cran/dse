\name{simple.monitoring}
\alias{simple.monitoring}
\alias{check.for.value.changes}
\alias{check.for.file.date.changes}
\alias{watch.data}
 
\title{Simple Monitoring}
\description{Automatic monitoring with e-mail of results.}
\usage{
	simple.monitoring(model, data.names, 
   previous.data=NULL,
   mail.list=NULL,
   error.mail.list=user.name(),
   message.title="Simple Monitoring",
   message.subject="Simple Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables= series.names(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   forecast.tag="f",
   run.again=F,
   save.as=NULL)
}
\arguments {
	\item{x}{see details below.}
}
\value{
	Invisibly return latest data for next comparison. This function is
	run mainly for its side effects.}
\details{
	model is a TSmodel. data.names is a TSdata (names) object.
	mail.list and error.mail.list should be single strings (not vectors)
	but the string can contain multiple user ids for mail
	If mail.list is NULL (default) then mail is not sent (useful for testing). 
	If error.mail.list  is  NULL   then mail is not sent (useful for testing).
	The default for error.mail.list is the result of user.name().
	This version does not allow for 
	    -combining forecasts (ie. monitoring.data or overriding data)
	    -input (policy) projections
	See combination.monitoring for these features.
	report.variables indicates output variables which are reported. It should
	correspond to a subset of names returned by series.names.
	data.tag and forecast.tag are tags placed beside data points in the report.
	If run.again is T then the monitoring is run regardless of data updates.
	show.end is the number of periods after the end of data (ie forecasts)
	 which should be displayed.
	show.start is the number of periods before the end of data (ie history)
	 which should be displayed. It is added to the end so it should be negative.
}
\sideffects{Mail is sent.}
\seealso{
	\code{\link{combination.monitoring}}
	\code{\link{simple.monitor.function.tests}}
}

\keyword{DSE}
\keyword{ts}


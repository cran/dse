\name{TSPADIdata2}
\alias{TSPADIdata2}
\title{Constructor for TSPADIdata object}
\description{
Make a TSPADIdata object from identifier information}
\usage{
	TSPADIdata2 <- function(input=NULL, output=NULL,
	  start = NA, end = NA, frequency = NA,
	  pad.start = FALSE, pad.end = FALSE,
	  start.server = NULL, server.process = NULL, cleanup.script = NULL,
	  stop.on.error = T, warn = T)
}
\arguments{
\item{ output  input}
{Lists of character vectors. Each vector has 5 elements indicating the 
server, the database, the series identifier, any tranformation, and a name.
At least one of input or output should be specified.}

\item{ start, end, frequency, pad.start, pad.end}
{ Information about data retrieval. See tfPADIdata for more details.}

\item{ start.server, server.process, cleanup.script, stop.on.error, warn}
{ Additional information server operation. See tfPADIdata for more details.}

}
\value{An object of class "TSPADIdata" "TSdata".}
\details{
This is an alternate way to produce an object which can be used for PADI database access.
Input and output values are passed to tfPADIdata so the TSdata object return is 
a list of tfPADIdata objects called input and output (but either of these can 
be NULL).
}
\seealso{
\code{\link{freeze.TSPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{tfPADIdata}}
\code{\link{TSPADIdata2}}
\code{\link{set.TSPADIdata}}
}
\examples{
inflation.sa.names <- TSPADIdata2(
	output = list(
         c("ets","", "ippia401",  "ytoypc", "industrial product price index"),
         c("ets","", "b820600",  "ytoypc", "cpi"),      
         c("ets","", "b820655",  "ytoypc", "cpixfe"),
         c("ets","", "b820672",  "ytoypc", "price goods excl food & energy"),
         c("ets","", "b820678",  "ytoypc", "price services"),
         c("ets","", "b800938",  "ytoypc", "total unit labour costs")
           ), stop.on.error=T, warn=T )

}
\keyword{DSE}
\keyword{ts}


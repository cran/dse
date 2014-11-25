\name{TSPADIdata}
\alias{TSPADIdata}
\alias{is.TSPADIdata}

\title{Constructor for TSPADIdata object}
\description{Make a TSPADIdata object from identifier information}
\usage{
	TSPADIdata( output=NULL,           input=NULL,
            output.transforms="",  input.transforms="", 
            output.names=NULL,     input.names=NULL,
            start=NA, end=NA, frequency=NA, 
            pad=FALSE, pad.start=pad, pad.end=pad,
            server="", db="",
	    start.server=NULL, server.process=NULL, cleanup.script=NULL,
            stop.on.error=T, warn=T)
  
	is.TSPADIdata(x)
}
\arguments{
\item{ output  input}
{Character vectors of series identifiers. At least one of input or output
should be specified.}

\item{ output.transforms input.transforms output.names input.names}
{ Character vectors of additional information corresponding to series 
identifiers. Scalar values are expanded to a vector of appropriate length.}

\item{ start, end, frequency, pad, pad.start, pad.end}
{ Information about data retrieval. See tfPADIdata for more details.}

\item{ server, db}
{ Information about the server. See tfPADIdata for more details.}

\item{ start.server, server.process, cleanup.script, stop.on.error, warn}
{ Additional information server operation. See tfPADIdata for more details.}
}
\value{An object of class "TSPADIdata" "TSdata".}
\details{
This is the way to produce an object which can be used for PADI database access.
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
	z <- TSPADIdata(output=c("seriesA", "seriesB"), server="myserver")
	is.TSPADIdata(z)
}
\keyword{DSE}
\keyword{ts}


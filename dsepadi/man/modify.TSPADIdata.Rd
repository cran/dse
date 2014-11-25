\name{modify.TSPADIdata}
\alias{modify.TSPADIdata}
\title{Modify a TSPADIdata Object}
\description{Modify a TSPADIdata object with new information}
\usage{
	modify.TSPADIdata(obj, output=NULL,  input=NULL,
            output.transforms="",  input.transforms="", 
            output.names=NULL,     input.names=NULL,
            start=NA, end=NA, frequency=NA, 
            pad.start=FALSE, pad.end=FALSE,
            server="", db="", start.server=NULL, 
            server.process=NULL, cleanup.script=NULL,
            stop.on.error=T, warn=T)
  
}

\arguments{
\item{ see TSPADIdata}{As for TSPADIdata.}
}
\value{An object of class "TSPADIdata" "TSdata".}
\details{This is the way to modify a TSPADIdata object.}
\seealso{
\code{\link{TSPADIdata}}
\code{\link{modify.tfPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{TSPADIdata}}
}
\examples{
inflation.sa.names <- TSPADIdata2(
	output = list(
         c("ets","", "b820678",  "ytoypc", "price services"),
         c("ets","", "b800938",  "ytoypc", "total unit labour costs")
           ), stop.on.error=T, warn=T )
z <- modify.TSPADIdata(inflation.sa.names, server="newserver")
}
\keyword{DSE}
\keyword{ts}


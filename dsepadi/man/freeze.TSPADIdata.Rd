\name{freeze.TSPADIdata}
\alias{freeze.TSPADIdata}
\title{Get fixed data snapshot using TSPADI}
\description{
Get a TSdata structure from TSPADI Database Interface}
\usage{
	freeze(data.id)
}
\arguments{
\item{ data.id}
{ An object of class TSPADIdata.}
}
\value{A TSdata object taken from the database.
}
\details{
This function extracts data using the TS PADI interface
to a database. See TSPADIdata, tfPADIdata, and freeze for more details. 
}
\seealso{
  
\code{\link{freeze}}
\code{\link{TSdata}}
\code{\link{TSPADIdata}}
}
\examples{
inflation.sa.names <- TSPADIdata2(
	output = list(
         c("ets","", "b820678",  "ytoypc", "price services"),
         c("ets","", "b800938",  "ytoypc", "total unit labour costs")
           ), stop.on.error=T, warn=T )
if(require("padi") & require("dse2")) {
z <- freeze(inflation.sa.names)}
}
\keyword{DSE}
\keyword{ts}


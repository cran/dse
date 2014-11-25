\name{retrieve.and.verify.data}
\alias{retrieve.and.verify.data}
\title{retrieve and verify data}
\description{Retrieve data from a database and verify a subset.}
\usage{
	retrieve.and.verify.data(data.names, 
	   verification.data=verification.data, fuzz=1e-10)
}
\arguments{
\item{ data.names}
{An object of class TSdata indicating databases and series.}

\item{ verification.data}{An object of class TSdata.}
\item{ fuzz}
{A (small) real number indicating the size of errors to tolerate.}
}
\value{An object of class TSdata containing data from the data bases is returned.
}
\sideffects{
If errors larger than fuzz are found then some warning messages are printed
and the series in question can be plotted.
}
\details{
The data from the data base is compared with verification.data, which would
normally be a subset (e.g. a subset used for model estimation). The
comparison is done  on the sub-sample for which verification data is
available. This provides
a mechanism for detecting large data revisions or re-basing.
}
\seealso{
\code{\link{TSPADIdata}}
\code{\link{set.TSPADIdata}}
\code{\link{set.tfPADIdata}}
\code{\link{freeze}}
}
\examples{
if(require("padi") & require("dse2")) {
   dnames <- TSPADIdata(output="D1", server="ets")
   d <- freeze(dnames)
   retrieve.and.verify.data(dnames, verification.data=d)}
}
\keyword{DSE}
\keyword{ts}


\name{ytoypc}
\alias{ytoypc}
\title{Convert to year to year percent change}
\description{Convert level data to year over year percent change. }
\usage{
    ytoypc(ser)}
\arguments{
    \item{ser}{A time series.}
}
\value{
A vector time series of the year over year percent change. This uses 
percent.change with lag=frequency(ser).
}
\seealso{
\code{\link{percent.change}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
z <- ytoypc(output.data(eg1.DSE.data))
}
%\keyword{DSE}
\keyword{ts}


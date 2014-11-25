\name{combine.TSdata}
\alias{combine.TSdata}
\title{Combine series from two TSdata objects.}
\description{Combine series from two TSdata objects.}
\usage{
	combine(e1, e2)}
\arguments{
	\item{e1, e2}{TSdata objects.}
}
\value{
An object of class TSdata which includes series from both e1 and e2.
}
\seealso{
tbind 
}
\examples{
if(is.R()) {data("eg1.DSE.data.diff", package="dse1")
            data("eg1.DSE.data", package="dse1") }
new.data.set <- combine(eg1.DSE.data.diff, eg1.DSE.data) 
}
\keyword{DSE}
\keyword{ts}


\name{percent.change}
\alias{percent.change}
\alias{percent.change.default}
\alias{percent.change.TSdata}
\alias{percent.change.list}
\alias{percent.change.TSestModel}
 
\title{Calculate percent change}
\description{Calculate the percent change relative to the data lag periods prior. 
}
\usage{
	percent.change(obj, base=NULL, lag=1, cumulate=F, e=F)}
\arguments{
\item{obj}
{An object of class TSdata or TSestModel, a time series matrix,
a matrix with columns corresponding to series (which are treated 
individually), or a list of one of these kinds of objects. (called m below)}
\item{e}
{If e is T the exponent of the series is used (after cumulating 
if cumulate is T). e can be
a logical vector with elements corresponding to columns of m.}
\item{ base}
{If base is provided it is treated as the first period value 
(that is, prior to differencing). It is prefixed to the m prior to 
cumulating. It should be a vector of length dim(m)[2]. 
(If e is T then base should be log of the original data).}
}
\value{
For an object of class TSdata the percent change calculation is done
with the output data and the result is an object of class TSdata (or a list of
objects of class TSdata).
For an object of class TSestModel the percent change calculation is done
with \code{estimates$pred} and the result is an object of class TSdata (or a 
list of objects of class TSdata).
}
\details{
If cumulate is T then the data is cumulated first. cumulate can be
a logical vector with elements corresponding to columns of m.
}
\seealso{
\code{\link{ytoypc}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
z <- percent.change(output.data(eg1.DSE.data))
}
\keyword{DSE}
\keyword{ts}


\name{percentChange}
\alias{percentChange}
\alias{percentChange.default}
\alias{percentChange.TSdata}
\alias{percentChange.TSestModel}
 
\title{Calculate percent change}
\description{Calculate the percent change relative to the data lag periods prior. 
}
\usage{
    percentChange(obj, ...)
    \method{percentChange}{default}(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE, ...)
    \method{percentChange}{TSdata}(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE, ...)
    \method{percentChange}{TSestModel}(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE, ...)
    }
\arguments{
    \item{obj}{
    An object of class TSdata or TSestModel, a time series matrix,
a matrix with columns corresponding to series (which are treated 
individually), or a list of one of these kinds of objects. (called m below)}
    \item{e}{
    If e is TRUE the exponent of the series is used (after cumulating 
if cumulate is TRUE). e can be
a logical vector with elements corresponding to columns of m.}
    \item{base}{
      If base is provided it is treated as the first period value 
      (that is, prior to differencing). It is prefixed to the m prior to 
      cumulating. It should be a vector of length dim(m)[2]. 
      (If e is TRUE then base should be log of the original data).}
    \item{lag}{integer indicating the number of periods relative to which the
      percent change should be calculated.}
    \item{cumulate}{logical indicating if the series should be cumulated before
      the percent change is calculated.}
    \item{...}{arguments passed to other methods.}
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
  Calculate the percent change relative to the data lag periods prior.
  obj should be  a  matrix or vector.
 
  If code{cumulate} is TRUE then the data is cumulated first. \code{cumulate} can be
  a logical vector with elements corresponding to columns of m.

  For a TSmodel the percent change calculation is done
  with input and output and the result is an object of class TSdata.

  For a TSestModel the percent change calculation is done
  with $estimates$pred and the result is an object of class TSdata
  }
\seealso{
\code{\link{ytoypc}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
z <- percentChange(outputData(eg1.DSE.data))
}
\concept{DSE}
\keyword{ts}


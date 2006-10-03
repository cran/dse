\name{percentChange}
\alias{percentChange}
\alias{percentChange.default}
 
\title{Calculate percent change}
\description{Calculate the percent change relative to the data lag periods prior. 
}
\usage{
    percentChange(obj, ...)
    \method{percentChange}{default}(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE, ...)
    }
\arguments{
    \item{obj}{
    An object on which the calculation is to be done. The default method works
    for a time series vector or matrix (with columns corresponding to series, which are treated 
    individually).}
    \item{e}{
    If e is TRUE the exponent of the series is used (after cumulating 
    if cumulate is TRUE). e can be
    a logical vector with elements corresponding to columns of obj.}
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
\value{The default method returns a time series matrix.}
\details{
  Calculate the percent change relative to the data lag periods prior.
  If \code{cumulate} is TRUE then the data is cumulated first. \code{cumulate} can be
  a logical vector with elements corresponding to columns of obj.
  }
\seealso{
\code{\link{ytoypc}}
}
}
\keyword{ts}

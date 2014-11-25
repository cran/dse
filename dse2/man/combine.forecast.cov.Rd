\name{combine.forecast.cov}
\alias{combine.forecast.cov}
\alias{combine.forecast.cov.estimators.wrt.data}
\alias{combine.forecast.cov.estimators.wrt.true} 

\title{Combine 2 Forecast Cov Objects}
\description{Combine 2 forecast.cov type objects.}
\usage{
    combine(obj1, obj2)
}
\arguments{
    \item{obj1, obj2}{
    Objects as returned by functions which calculate forecast covariances.}
}
\value{
An object as returned by functions which calculate forecast covariances.
}
\details{
Functions which calculate forecast covariances return lists. Usually multiple 
estimation techniques or models will be combined together when the object is
first formed. However, it is sometimes useful to add results calculated later
without re-doing the initial object.
}
\seealso{
\code{\link{forecast.cov.estimators.wrt.data}},
\code{\link{forecast.cov.estimators.wrt.true}}
\code{\link{forecast.cov}}
}
\examples{#z <-  combine(obj1, obj2)}
%\keyword{DSE}
\keyword{ts}


\name{trimNA}
\alias{trimNA}
\alias{trimNA.default}

\title{Trim NAs from Time Series}
\description{Trim NAs from the start and end of a time series object.}
\usage{
    trimNA(x, start. =TRUE, end. =TRUE)
    \method{trimNA}{default}(x, start. =TRUE, end. =TRUE)
    }
\arguments{
    \item{x}{A time series matrix or an object of class TSdata.}
    \item{start.}{If start.=FALSE then beginning NAs are not trimmed.}
    \item{end.}{If end.=FALSE   then ending NAs are not trimmed.}
}
\value{A time series object which is windowed to the subset of 
 data which is not NAs (usually the available data).}
\details{
Trim NAs from the ends of a time series object.
Observations in a given period for all series are dropped if any one contains an NA.
}
\examples{
trimNA(ts(rbind(NA, matrix(1:20,10,2)), start=c(1980,1), frequency=12))
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



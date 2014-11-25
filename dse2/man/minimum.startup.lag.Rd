\name{minimum.startup.lag}
\alias{minimum.startup.lag}
\alias{minimum.startup.lag.SS}
\alias{minimum.startup.lag.ARMA}
\alias{minimum.startup.lag.TSestModel}
\alias{start.shift}

\title{Starting Periods Required}
\description{Number of Starting Periods Required for a Model}
\usage{
    minimum.startup.lag(model)
}
\arguments{
    \item{model}{A TSmodel or object containing a TSmodel.}
}
\value{An integer.}
\details{
    For many time series models several starting data points are 
    required before output from the model can be calculated (or makes
    sense). This generic function extracts or calculates that 
    number of periods.
}

\seealso{
    \code{\link{TSmodel}}
}

\note{There is some redundancy between this and start.shift which 
      should be cleaned up.}
%\keyword{DSE}
\keyword{ts}


\name{input.dimension}
\alias{input.dimension}
\alias{input.dimension.default}
\alias{input.dimension.TSdata}
\alias{input.dimension.SS}
\alias{input.dimension.ARMA}
\alias{input.dimension.TSestModel}

\alias{output.dimension}
\alias{output.dimension.default}
\alias{output.dimension.TSdata}
\alias{output.dimension.SS}
\alias{output.dimension.ARMA}
\alias{output.dimension.TSestModel}

\title{Number of Series in in Input or Output}
\description{Number of input or output series in a TSdata object.}
\usage{
    input.dimension(x)
    \method{input.dimension}{default}(x)
    \method{input.dimension}{TSdata}(x)
    \method{input.dimension}{SS}(x)
    \method{input.dimension}{ARMA}(x)
    \method{input.dimension}{TSestModel}(x)
    
    output.dimension(x)
    \method{output.dimension}{default}(x)
    \method{output.dimension}{TSdata}(x)
    \method{output.dimension}{SS}(x)
    \method{output.dimension}{ARMA}(x)
    \method{output.dimension}{TSestModel}(x)
}
\arguments{
    \item{x}{Object of class TSdata, TSmodel or TSestModel.}
}
\value{
    An integer indicating the number of series.
}
\seealso{
    \code{\link{input.series.names}}
    \code{\link{output.series.names}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
output.dimension(eg1.DSE.data) 
}
%\keyword{DSE}
\keyword{ts}


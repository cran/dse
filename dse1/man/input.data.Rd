\name{input.data}
\alias{input.data}
\alias{input.data.default}
\alias{input.data.TSdata}
\alias{input.data.TSestModel}
\alias{input.data<-}
\alias{input.data<-.default}
\alias{input.data<-.TSdata}
\alias{input.data<-.TSestModel}

\alias{output.data}
\alias{output.data.default}
\alias{output.data.TSdata}
\alias{output.data.TSestModel}
\alias{output.data<-}
\alias{output.data<-.default}
\alias{output.data<-.TSdata}
\alias{output.data<-.TSestModel}

\title{TSdata Series}
\description{Extract or set input or output series in a TSdata object.}
\usage{
    input.data(x, series=seqN(input.dimension(x)))
    \method{input.data}{default}(x, series=seqN(input.dimension(x)))
    \method{input.data}{TSdata}(x,  series=seqN(input.dimension(x)))
    \method{input.data}{TSestModel}(x, series=seqN(input.dimension(x)))

    output.data(x, series=seqN(output.dimension(x)))
    \method{output.data}{default}(x, series=seqN(output.dimension(x)))
    \method{output.data}{TSdata}(x,  series=seqN(output.dimension(x)))
    \method{output.data}{TSestModel}(x, series=seqN(output.dimension(x)))

    input.data(data)  <- x
    output.data(data) <- y
}
\arguments{
    \item{data}{object of class TSdata.}
    \item{x}{a time series matrix.}
    \item{y}{a time series matrix.}
    \item{series}{vector of strings or integers indicating the series to select.}
}
\value{
The first usages returns the input or output series.
The second usages assigns the input or output series.
}
\seealso{
    \code{\link{TSdata}}
    \code{\link{select.series}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
output.data(eg1.DSE.data) 
}
%\keyword{DSE}
\keyword{ts}


\name{input.series.names}
\alias{input.series.names}
\alias{input.series.names.TSdata}
\alias{input.series.names.TSmodel}
\alias{input.series.names.TSestModel}
\alias{input.series.names<-}
\alias{input.series.names<-.TSdata}
\alias{input.series.names<-.TSmodel}
\alias{input.series.names<-.TSestModel}


\alias{output.series.names}
\alias{output.series.names.TSdata}
\alias{output.series.names.TSmodel}
\alias{output.series.names.TSestModel}
\alias{output.series.names<-}
\alias{output.series.names<-.TSdata}
\alias{output.series.names<-.TSmodel}
\alias{output.series.names<-.TSestModel}

\title{TSdata Series Names}
\description{Extract or set names of input or output series in a TSdata object.}
\usage{
    input.series.names(x)
    \method{input.series.names}{TSdata}(x)
    \method{input.series.names}{TSmodel}(x)
    \method{input.series.names}{TSestModel}(x)

    output.series.names(x)
    \method{output.series.names}{TSdata}(x)
    \method{output.series.names}{TSmodel}(x)
    \method{output.series.names}{TSestModel}(x)

    input.series.names(x)  <- value
    output.series.names(x) <- value
}
\arguments{
    \item{x}{Object of class TSdata, TSmodel or TSestModel.}
}
\value{
The first usages gives a vector of strings with the series names.
The second usages assigns a vector of strings to be the series names of data.
}
\seealso{
    \code{\link{series.names}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
output.series.names(eg1.DSE.data) 
}
%\keyword{DSE}
\keyword{ts}



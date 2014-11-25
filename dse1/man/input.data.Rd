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
	input.data(data)
	output.data(data)
	input.data(data)  <- x
	output.data(data) <- y
}
\arguments{
	\item{data}{Object of class TSdata.}
	\item{x}{A time series matrix.}
	\item{y}{A time series matrix.}
}
\value{
The first usages returns the input or output series.
The second usages assigns the input or output series.
}
\seealso{
	\code{\link{TSdata}}
}
\examples{
if(is.R()) data("eg1.DSE.data", package="dse1")
output.data(eg1.DSE.data) 
}
\keyword{DSE}
\keyword{ts}


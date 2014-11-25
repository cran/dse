\name{input.periods}
\alias{ input.periods}
\alias{output.periods}
\alias{ input.periods.TSdata}
\alias{output.periods.TSdata}
\alias{ input.periods.TSestModel}
\alias{output.periods.TSestModel}

\alias{ input.start}
\alias{output.start}
\alias{ input.start.TSdata}
\alias{output.start.TSdata}

\alias{ input.end}
\alias{output.end}
\alias{ input.end.TSdata}
\alias{output.end.TSdata}

\alias{ input.frequency}
\alias{output.frequency}
\alias{ input.frequency.TSdata}
\alias{output.frequency.TSdata}

\title{TSdata Periods}
\description{Apply a method to the input or output data.}
\usage{
    input.periods(x, ...)
    \method{input.periods}{TSdata}(x)
    \method{input.periods}{TSestModel}(x)

    output.periods(x, ...)
    \method{output.periods}{TSdata}(x)
    \method{output.periods}{TSestModel}(x)

    input.start(x)
    \method{input.start}{TSdata}(x)
    \method{input.start}{TSestModel}(x)

    output.start(x)
    \method{output.start}{TSdata}(x)
    \method{output.start}{TSestModel}(x)

    input.end(x)
    \method{input.end}{TSdata}(x)
    \method{input.end}{TSestModel}(x)

    output.end(x)
    \method{output.end}{TSdata}(x)
    \method{output.end}{TSestModel}(x)

    input.frequency(x)
    \method{input.frequency}{TSdata}(x)
    \method{input.frequency}{TSestModel}(x)

    output.frequency(x)
    \method{output.frequency}{TSdata}(x)
    \method{output.frequency}{TSestModel}(x)

}
\arguments {
    \item{x}{An object containing TSdata.}
}
\value{Depends.}
\details{
    Apply a method to the input or output data so, for example,
    input.periods(x) in theory does periods(input.data(x)), which 
    returns the number of periods in input data. The actual implementation
    may not do periods(input.data(x)). For example, with TSPADIdata 
    input.data(x) requires a database retrieval which may be fairly slow,
    while the number of periods may be available much more quickly.
}

\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    output.periods(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


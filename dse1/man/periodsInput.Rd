\name{periodsInput}
\alias{periodsInput}
\alias{periodsInput.TSdata}
\alias{periodsInput.TSestModel}

\alias{periodsOutput}
\alias{periodsOutput.TSdata}
\alias{periodsOutput.TSestModel}

\alias{input.start}
\alias{input.start.TSdata}
\alias{input.start.TSestModel}

\alias{output.start}
\alias{output.start.TSdata}
\alias{output.start.TSestModel}

\alias{input.end}
\alias{input.end.TSdata}
\alias{input.end.TSestModel}

\alias{output.end}
\alias{output.end.TSdata}
\alias{output.end.TSestModel}

\alias{input.frequency}
\alias{input.frequency.TSdata}
\alias{input.frequency.TSestModel}

\alias{output.frequency}
\alias{output.frequency.TSdata}
\alias{output.frequency.TSestModel}

\title{TSdata Periods}
\description{Apply a method to the input or output data.}
\usage{
    periodsInput(x)
    \method{periodsInput}{TSdata}(x)
    \method{periodsInput}{TSestModel}(x)

    periodsOutput(x)
    \method{periodsOutput}{TSdata}(x)
    \method{periodsOutput}{TSestModel}(x)

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
\arguments{
    \item{x}{An object containing TSdata.}
}
\value{Depends.}
\details{
    Apply a method to the input or output data so, for example,
    periodsInput(x) in theory does periods(input.data(x)), which 
    returns the number of periods in input data. The actual implementation
    may not do periods(input.data(x)). For example, with TSPADIdata 
    input.data(x) requires a database retrieval which may be fairly slow,
    while the number of periods may be available much more quickly.
}

\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    periodsOutput(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


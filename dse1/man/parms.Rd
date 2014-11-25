\name{parms}
\alias{parms}
\alias{parms.TSmodel}
\alias{parms.TSestModel}

\title{Extract Model Parameters.}
\description{
    A generic function to extract parameter values from some DSE objects.
}
\usage{
    parms(model)
    \method{parms}{TSmodel}(model)
    \method{parms}{TSestModel}(model)
    }
\arguments{    \item{obj}{An object of class TSmodel or TSestModel.}}
\value{A vector of parameter values.}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    parms(model)
}
%\keyword{DSE}
\keyword{ts}


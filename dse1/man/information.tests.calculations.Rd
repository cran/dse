\name{information.tests.calculations}
\alias{information.tests.calculations}

\title{Calculate selection criteria}
\description{Calculates several model selection criteria.}
\usage{
    information.tests.calculations(lst, sample.start=1,sample.end=NULL, warn=T)
}
\arguments{
    \item{lst}{One or more objects of class TSestModel.}
    \item{sample.start}{
        The start of the period to use for criteria calculations.
    }
    \item{sample.end}{
        The end of the period to use for criteria calculations. 
        If omitted the end of the sample is used.
    }
    \item{warn}{If FALSE then some warning messages are suppressed.}
}
\value{The calculated values are returned in a vector.}
\seealso{
    \code{\link{information.tests}}
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    information.tests.calculations(model)
}
%\keyword{DSE}
\keyword{ts}


\name{McMillan.degree}
\alias{McMillan.degree}
\alias{McMillan.degree.ARMA}
\alias{McMillan.degree.SS}
\alias{McMillan.degree.TSestModel}

\title{Calculate McMillan Degree}
\description{

Calculate the McMillan degree of an ARMA TSmodel.
}
\usage{
    McMillan.degree(model,  ...)
    \method{McMillan.degree}{ARMA}(model, fuzz=1e-4, verbose=TRUE, warn=TRUE, ...)
    \method{McMillan.degree}{SS}(model, fuzz=1e-4, ...)
    \method{McMillan.degree}{TSestModel}(model,  ...)
}
\arguments{
    \item{model}{An object of class TSmodel.}
    \item{fuzz}{Roots within fuzz distance are counted as equivalent.}
    \item{verbose}{If TRUE roots are printed.}
    \item{warn}{If FALSE then warnings about unit roots added for TREND are not printed.}
    \item{...}{arguments to be passed to other methods.}
}
\value{
A list with elements gross and distinct containing all roots 
and distinct roots.
}
\sideffects{
The number of roots and distinct roots is printed if verbose is TRUE.
}
\seealso{

\code{\link{stability}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
McMillan.degree(model)
}
%\keyword{DSE}
\keyword{ts}


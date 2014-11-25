\name{roots}
\alias{roots}
\alias{roots.SS}
\alias{roots.ARMA}
\alias{roots.TSestModel}

\title{Calculate Model Roots}
\description{Calculate roots of a TSmodel.}
\usage{
    roots(obj, ...)
    \method{roots}{SS}(obj, fuzz=0, randomize=FALSE, ...)
    \method{roots}{ARMA}(obj, fuzz=0, randomize=FALSE, warn=TRUE, by.poly=FALSE, ...)
    \method{roots}{TSestModel}(obj, ...)
    }
\arguments{
    \item{obj}{An object of class TSmodel.}
    \item{fuzz}{If non-zero then roots within fuzz 
       distance are considered equal.}
    \item{randomize}{
       Randomly arrange complex pairs of roots so the one with the positive imaginary
       part is not always first (so random experiments are not biased).}
    \item{warn}{If FALSE then warnings about unit roots added for TREND are not printed.}
    \item{by.poly}{
       If TRUE then roots are calculated by expanding the determinant of the A
       polynomial.  Otherwise, they are calculated by converting to a state
       space representation and calculating the eigenvalues of F. This second
       method is preferable for speed, accuracy, and because of a limitation
       in the degree of a polynomial which can be handled by polyroot.}
    \item{verbose}{arguments to be passed to other methods.}
    \item{...}{arguments passed to other methods.}
}
\value{

The eigenvalues of the state transition matrix or the inverse of the roots of the
determinant of the AR polynomial are returned.
}
\seealso{
\code{\link{stability}},
\code{\link{McMillan.degree}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- est.VARX.ls(eg1.DSE.data.diff)
roots(model)
}
%\keyword{DSE}
\keyword{ts}


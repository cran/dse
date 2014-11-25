\name{findg}
\alias{findg}
\title{Find Equivalence Transformation}
\description{

Try to find a matrix g which makes two models equivalent.
}
\usage{
    findg(model1, model2, minf=nlmin)
}
\arguments{

    \item{model1, model2}{Objects of class TSmodel.}
    \item{minf}{Algorithm used to minimize the difference.}
}
\value{A matrix which converts one model to the other.}
\details{
This is set up as a minimization problem with the objective 
to reduce the squared difference between parameters.
}
\note{WARNING: This program does not work very well.}
\seealso{
    \code{\link{gmap}}
}
\examples{
# findg(model1, model2)
}
%\keyword{DSE}
\keyword{ts}


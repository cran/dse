\name{generate.model.SS}
\alias{generate.model.SS}
\title{Randomly generate a state space model}
\description{Randomly generate a state space model.}
\usage{
    generate.model.SS(m,n,p, stable=FALSE)
}
\arguments{
    \item{n,m,p}{Input, state and output dimensions.}
    \item{stable}{TRUE or FALSE indicating if the model must be stable.}
}
\value{An SS TSmodel.}
\details{
Randomly generate a state space model. If stable is true then the
largest root will have magnitude less than 1.0.
}
\examples{
z <-  generate.model.SS(2,3,1)
}
%\keyword{DSE}
\keyword{ts}


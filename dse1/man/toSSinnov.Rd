\name{toSSinnov}
\alias{toSSinnov}
\title{Convert to State Space Innovations Model}
\description{Convert to a state space innovations representation.}
\usage{
    toSSinnov(model, ...)}
\arguments{
    \item{model}{an object of class TSmodel.}
    \item{...}{arguments passed to other methods.}
}
\value{
An object of class 'SS' 'TSmodel' containing an 
innovations state space representation.
}
\seealso{
\code{\link{toSS}},
\code{\link{toSSOform}}
\code{\link{toSSChol}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- estVARXls(eg1.DSE.data.diff)
model <- toSSinnov(model)
}
%\keyword{DSE}
\keyword{ts}


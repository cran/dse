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
If the argument is a TSmodel then the result is an object of 
class 'SS' 'TSmodel' If the argument is a TSestModel then the converted model is
evaluated with the data an a TSestModel is returned. The TSmodel is an 
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
summary(model)

model2 <- SS(F=diag(1,3), H=matrix(c(1,0,0,1,0,0),2,3),
   Q=diag(0.5, 3, 3), R=diag(1.1, 2,2),
    description="test model", output.names=c("output 1", "output 2"))
model2 <- toSSinnov(model2)
summary(model2)
}
%\keyword{DSE}
\keyword{ts}


\name{tags}
\alias{tags}
\alias{tags<-}
\alias{tagged}
\alias{tagged.default}
\alias{tagged.TSdata}
\alias{is.tagged}

\title{Tagged Matrices}
\description{Matrices with an decription for each point. }
\usage{
	tags(x)
	tags(x) <- value
	tagged(x, value)
	tagged(x, input.value, output.value)
	is.tagged(obj)
}
\arguments {
	\item{x}{A matrix or TSdata object.
	}
	\item{value}{
		If x is a matrix then value should be a matrix of strings of 
		the corresponding dimension or a single string which is expanded
		to the dimension of x.
	\item{input.value}{
		If x is TSdata then input.value should be a matrix of strings of 
		the same dimension as input.data(x)or a single string which is 
		expanded.
	\item{output.value}{
		If x is TSdata then output.value should be a matrix of strings of 
		the same dimension as input.data(x)or a single string which is 
		expanded.
	}
}
\value{
	A matrix of class "tagged" or a TSdata object with specific 
	class "tagged"}
\details{
	The matrix of descriptive information is assigned as an attribute of the
	matrix. Most testing of the methods for this class has been with a 
	single character tag which can be used as a flag, for example, to 
	indicate the different sources for the data points.
}
\seealso{
	\code{\link{TSdata}}
}
\examples{
	x <- tagged(matrix(rnorm(100), 50,2), "r")
	is.tagged(x)
}

\keyword{DSE}
\keyword{ts}


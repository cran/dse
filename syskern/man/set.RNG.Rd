\name{set.RNG}
\alias{set.RNG}
\alias{RNGkind}
\alias{random.number.test}
\title{Set the random number generator}
\description{Set the random number generator.}
\usage{
	set.RNG(kind=NULL, seed=NULL, normal.kind=NULL)}
\arguments{None required} {
	\item{kind}{
		A character string.
	}
	\item{seed}{
		A vector of numbers (depending on kind).
	}
	\item{normal.kind}{
		A character string.
	}
}
\value{The old setting.}
\details{
	Sets the uniform and normal random number generators and the seed.
	The function \code{random.number.test()) can be used to check that the
	function is working properly. The old setting is returned in a format
	which can be used in another call to set.RNG. (This would reset to the
	original value.) If no arguments are given the current setting is 
	returned. In R see RNGkind for more details.
}
\sideffects{Sets global variables controlling the uniform and normal random 
	number generators and the global seed.}
\seealso{
	\code{\link{set.seed}}
	\code{\link{RNGkind}}
	\code{\link{runif}}
	\code{\link{rnorm}}
}
\examples{
	set.RNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
	rnorm(10)
	random.number.test()
}
\keyword{syskern}
\keyword{OS}
\keyword{S/R}


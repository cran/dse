\name{juice.function.tests}
\alias{juice.function.tests}
\alias{juice.graphics.tests}
\title{Test functions in Juice Package of DSE Bundle}
\description{Test proper functioning of part of the DSE library.}
\usage{
	juice.function.tests(verbose=T, synopsis=T, pause=F,
      		fuzz.small=1e-14, fuzz.large=1e-6)
}

\arguments {
	\item{verbose}{Print extra detail.}
	\item{synopsis}{Indicate if all tests passed or not.}
	\item{fuzz.small}{Tight tolerance numerical comparisons.}
	\item{fuzz.large}{Loose tolerance numerical comparisons.}
	\item{pause}{T or F indicating if the device should be 
		paused after each page of graphics}
}

\details{
	This function runs several tests. Some numerical tests typically can
	be checked at a tighter tolerance than others. Many of the tests rely
	on reproducible random number generation. (See set.RNG.)
}
\sideffects{Test results are printed and in some cases plots are produced.}
\examples{
	\dontrun{juice.function.tests()}
	\dontrun{juice.graphics.tests(pause=F)}
}
\keyword{DSE}
\keyword{time series}


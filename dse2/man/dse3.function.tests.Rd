\name{dse3.function.tests}
\alias{dse3.function.tests}
\alias{dse3.graphics.tests}
\title{Function Tests for Part of DSE Bundle}
\description{Test proper functioning of part of the DSE library.}
\usage{
	dse3.function.tests(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	dse3.graphics.tests(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
}

\arguments {
	\item{verbose}{Print extra detail.}
	\item{synopsis}{Indicate if all tests passed or not.}
	\item{fuzz.small}{Tight tolerance numerical comparisons.}
	\item{fuzz.large}{Loose tolerance numerical comparisons.}
}

\details{
	This function runs several tests. Some numerical tests typically can
	be checked at a tighter tolerance than others. Many of the tests rely
	on reproducible random number generation. (See set.RNG.)
}
\sideffects{Test results are printed and in some cases plots are produced.}
\examples{
	\dontrun{dse3.function.tests()}
	\dontrun{dse3.graphics.tests()}
}
\keyword{DSE}
\keyword{ts}


\name{dse1.function.tests}
\alias{dse1.function.tests}
\title{dse1.function.tests}
\description{Test proper functioning of part of the DSE library.}
\usage{
	dse1.function.tests(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)}

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
	dse1.function.tests()
}
\keyword{DSE}
\keyword{ts}


\name{dsecurvature.function.testsA}
\alias{dsecurvature.function.testsA}
\alias{dsecurvature.function.testsB}
\alias{dsecurvature.function.testsC}
\title{dsecurvature.function.tests}
\description{Test proper functioning of part of the DSE library.}
\usage{
	dsecurvature.function.testsA(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	dsecurvature.function.testsB(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	dsecurvature.function.testsC(verbose=T, synopsis=T,
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
	Don't attempt testsB unless you have a fast computer with lots of 
	memory. If you complete testsC please let me know (PG).
}
\sideffects{Test results are printed and in some cases plots are produced.}
\examples{
	\dontrun{dsecurvature.function.testsA()}
	\dontrun{dsecurvature.function.testsB()}
	\dontrun{dsecurvature.function.testsC()}
}
\keyword{DSE}
\keyword{ts}


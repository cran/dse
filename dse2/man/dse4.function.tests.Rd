\name{dse4.function.tests}
\alias{dse4.function.tests}
\alias{dse4.graphics.tests}
\title{Test Functions of Part of DSE Package}
\description{Test proper functioning of part of the DSE library.}
\usage{
	dse4.function.tests(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	dse4.graphics.tests(verbose=T, synopsis=T,
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
	\dontrun{dse4.function.tests(graphics=F)}
	\dontrun{dse4.graphics.tests(pause=F)}
}
\keyword{DSE}
\keyword{ts}


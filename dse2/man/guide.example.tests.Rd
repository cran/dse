\name{guide.example.tests}
\alias{guide.example.tests.part1}
\alias{guide.example.tests.part2}
\alias{guide.example.tests.part3}
\title{Users Guide Tests}
\description{Test examples in the DSE Users' Guide.}
\usage{
	guide.example.tests.part1(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	guide.example.tests.part2(verbose=T, synopsis=T,
      		fuzz.small=1e-14, fuzz.large=1e-6)
	guide.example.tests.part3(verbose=T, synopsis=T,
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
	\dontrun{guide.example.tests.part1()}
	\dontrun{guide.example.tests.part2()}
	\dontrun{guide.example.tests.part3()}
}
\keyword{DSE}
\keyword{ts}


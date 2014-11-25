\name{DSECOMPILED}
\alias{.DSECOMPILED}
\title{Use Compiled or S/R Code}
\description{
	Determines if compiled code should be used or not.
	}
\usage{
	.DSECOMPILED <- T
	.DSECOMPILED <- F
}

\sideffects{
	Several DSE functions which call compiled fortran or C code will use
	the equivalent S/R version if this global variable is set to FALSE.
	Its main purpose is for debugging.
}
\examples{
	.DSECOMPILED <- T
}
\keyword{DSE}
\keyword{ts}


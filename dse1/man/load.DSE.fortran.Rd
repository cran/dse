\name{load.DSE.fortran}
\alias{load.DSE.fortran}
\alias{dse.win.for.tab}
\title{Load compiled routines for time series models}
\description{Load compiled routines for time series models (Splus only).}
\usage{
	load.DSE.fortran(from=NULL, large=F)}
\arguments{
	\item{from}{Directory where the compiled object code is stored.}
	\item{large}{
	If very large states (>100) are required it is necessary to load 
	a different version of the compiled code, which requires compiling a second
	version. This makes the loaded code, and thus S, very large and so it should
	only be used if required.
	}
}
\value{invisible list of any errors}
\details{
	The default assumes the compiled object code is in a subdirectory
	named DSE.HOME/lib where DSE.HOME is a global variable set 
	in .First.lib or when DSE is installed. 
	This function is not needed in R.
}
\sideffects{compiled routines are loaded.}
\examples{if(is.S()) load.DSE.fortran()}
\keyword{DSE}
\keyword{ts}


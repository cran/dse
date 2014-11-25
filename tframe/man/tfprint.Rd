\name{tfprint}
\alias{tfprint}
\alias{tfplot}
\alias{tfplot.default}
\alias{tfOnePlot}
\alias{tfprint.default}
\alias{tfprint.tframe}
\alias{tfprint.tframe.default}

\title{Print and Plot Tframed Objects}
\description{
	Print or plot tframe or tframed objects. 
}
\usage{
	tfprint(x)
	tfplot(x)
}

\arguments{
	\item{x}{A tframe or tframed object.}
}
\value{tfprint methods return the object invisibly.}
\details{
	In many cases these are the 
	same as print and plot methods. However, those methods are defined 
	in different libraries. Also, new classes of time series may define
	there own print and plot methods in ways which use a different logic
	from the tframe library. Thus tfprint and tfplot provide another
	way to program functions which use methods consistent with the tframe
	library logic.
	
	Note that the start and end arguments to tfplot are "start." and "end.".
	In some situations truncation works and these can be abreviated, but
	often it does not work and the error message is opague.
}
\sideffects{An object is printed or plotted.}
\seealso{
	\code{\link{tframe}}
	\code{\link{tframed}}
	\code{\link{print}}
	\code{\link{plot}}
}
\examples{
	tfplot(ts(rnorm(100), start=c(1982,1), frequency=12))
	tfplot(ts(rnorm(100), start=c(1982,1), frequency=12), start.=c(1985,6))
}
\keyword{tframe}
\keyword{ts}


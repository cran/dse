\name{is.unix}
\alias{is.unix}
\alias{is.Linux}
\alias{is.Sun4}
\alias{is.Sun5}
\alias{is.Win3.1}
\alias{is.Win32}
\alias{is.Win95}
\alias{is.WinNT}
\alias{is.Mac}
\alias{is.MSwindows}

\title{Determine the Operating System}
\description{Determine the operating system.}

\usage{
	is.unix()
	is.Linux()
	is.Sun4()
	is.Sun5()
	is.Win3.1()
	is.Win32()
	is.Win95()
	is.WinNT()
	is.Mac()
	is.MSwindows()
}

\arguments{None}
\value{TRUE or FALSE}
\details{
	If the operating system is the one in question then TRUE is
	returned, otherwise FALSE. Windows and Mac versions are not yet
	done. (Offers welcome.)
}
\seealso{
	\code{\link{is.R}}
	\code{\link{is.S}}
}
\examples{is.unix() }
\keyword{syskern}
\keyword{OS}
\keyword{S/R}


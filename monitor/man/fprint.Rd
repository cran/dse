\name{fprint}
\alias{fprint}
\alias{fprint.tagged}
\title{Formatted Printing of Time Series}
\description{Generate a formatted character matrix of time series data.}
\usage{
	fprint(matrix, super.title=NULL, sub.title=NULL,
          digits=options()$digits, space=" ", file=NULL, append=F)
}
\arguments {
	\item{matrix}{A matrix of time series data.}
	\item{super.title}{Top heading.}
	\item{sub.title}{Second level of headings.}
	\item{digits}{Number of digits to print.}
	\item{file}{Name of a file to generate.}
	\item{append}{If T output is appended to the file, otherwise 
	  the file is overwritten.}
}
\value{A character matrix.}
\details{
	This is a generic method for print a table of time series data with
	formatting control. Currently the only specific method is fprint.tagged. 
}
\sideffects{If \code{file} is specified then a file will be created.}
\seealso{
	\code{\link{print}}
	\code{\link{tfprint}}
}
\keyword{DSE}
\keyword{ts}
\keyword{programming}
\keyword{utilities}


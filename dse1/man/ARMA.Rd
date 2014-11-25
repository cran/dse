\name{ARMA}
\alias{ARMA}
\alias{is.ARMA}

\title{ARMA Model Constructor}
\description{Constructs an ARMA TSmodel object as used by the DSE package.}
\usage{
	ARMA(A=NULL, B=NULL, C=NULL, TREND=NULL, description=NULL,
          names=NULL, input.names=NULL, output.names=NULL)
	is.ARMA(...)
}
\arguments {
	\item{A}{(axpxp) array. (The auto-regressive polynomial array.)}
	\item{B}{(bxpxp) array. (The moving-average polynomial array.)}
	\item{C}{(cxpxm) array. (The  input polynomial array.) C should be NULL if 
	there is no input}
}
\value{An ARMA TSmodel}
\details{ 
The ARMA model is defined by:

      Ay =  Bw + Cu
      
sometimes written

      A(L)y(t) =  B(L)w(t) + C(L)u(t)

where
	\item{A}{(axpxp) is the auto-regressive polynomial array.}
	\item{B}{(bxpxp) is the moving-average polynomial array.}
	\item{C}{
	(cxpxm) is the  input polynomial array. C should be NULL if 
	there is no input}
	\item{y}{is the p dimensional output data.}
	\item{u}{is the m dimensional control (input) data.}

}

\seealso{
	\code{\link{TSmodel}}
}
\examples{
	mod1 <- ARMA(A=array(c(1,-.25,-.05), c(3,1,1)), B=array(1,c(1,1,1)))
}

\keyword{DSE}
\keyword{ts}


\name{acf.M}
\alias{acf.M}
\alias{acf.M.TSdata}
\alias{acf.M.TSmodel}
\alias{acf.M.TSestModel}

\title{Calculate auto-covariance}
\description{

Calculate a matrix with partitions [M0|...|Mi|...|Ml] giving the auto-covariance.
}
\usage{
	acf.M(object, type ='covariance', sub.mean=T, Psi=NULL, lag=)
}
\arguments{

\item{ object}
{
An object of class TSdata or TSmodel.}

}


\item{ type}
{
With the defaults the blocks are auto-covariances. If type == 'correlation' 
the result is scaled to give autocorrelations.
\item{ sub.means}
{
Only valid if object is of class TSdata. If F then means are not subtracted.
\item{ Psi}
{
A matrix of innovation covariance. Only valid if object is of class TSmodel.

}
\value{

A matrix with partitions [M0|...|Mi|...|Ml] giving the covariance or correlation, 
including the that between the output and input series 
(as in the first block row of a Hankel matrix).
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- acf.M(eg1.DSE.data.diff)
model <- TSmodel(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
#  z <- acf.M(model) not working
}
\keyword{DSE}
\keyword{ts}


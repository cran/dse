\name{acf.M}
\alias{acf.M}
\alias{acf.M.TSdata}
\alias{acf.M.TSmodel}
\alias{acf.M.TSestModel}

\title{Calculate Auto-Covariance}
\description{
    Calculate a matrix with partitions [M0|...|Mi|...|Ml] giving 
    the auto-covariance.
}
\usage{
    acf.M(obj, ...)
    \method{acf.M}{TSdata}(obj, lag=round(6*log(periods(obj))), 
        type ="covariance", sub.mean=T)
    \method{acf.M}{TSmodel}(obj, lag=NULL, type ="covariance", Psi=NULL)
    \method{acf.M}{TSestModel}(obj, ...)
}
\arguments{
    \item{obj}{An object of class TSdata or TSmodel.}
    \item{lag}{Number of lags for which to calculate the autocorrelations.}
    \item{type}{
        With the defaults the blocks are auto-covariances. 
        If type == 'correlation' the result is scaled to give 
        autocorrelations.}
    \item{sub.mean}{
        Only valid if object is of class TSdata. If F then means 
        are not subtracted.}
    \item{Psi}{
        A matrix of innovation covariance. Only valid if object 
        is of class TSmodel.}
    \item{...}{arguments to be passed to other methods.}
}

\value{
    A matrix with partitions [M0|...|Mi|...|Ml] giving the covariance 
    or correlation, including the that between the output and input 
    series (as in the first block row of a Hankel matrix).
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    z <- acf.M(eg1.DSE.data.diff)
    model <- TSmodel(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
    #  z <- acf.M(model) not working
}
%\keyword{DSE}
\keyword{ts}


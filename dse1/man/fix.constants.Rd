\name{fix.constants}
\alias{fix.constants}
\title{Fix TSmodel Coefficients (Parameters) to Constants}
\description{
    Fix any coefficients within fuzz of 0.0 or 1.0 to exactly 0.0 or 1.0.
    This will not change the model much but will affect some estimation 
    techniques and information criteria results, as these are considered to 
    be constants rather than coefficients.
}
\usage{
    fix.constants(model, fuzz=1e-5, constants=NULL)
}
\arguments{
    \item{model}{an object of class TSmodel.}
    \item{fuzz}{absolute difference to be considered equivalent.}
    \item{constants}{NULL or a list of logical arrays.}
}
\value{
    An object of class 'SS' 'TSmodel' with parameters near 0.0 or 1.0 set 
    as constants 0.0 or 1.0. If constants is not NULL then it should be a list 
       with logical (T/F) arrays named const.F, const.G ..., for any 
       arrays in which there are elements (not == 0 or 1) which are to 
       be treated as constant.
}

\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- to.ARMA(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
    model <- fix.constants(model)
}
\seealso{
    \code{\link{fixF}}
}
%\keyword{DSE}
\keyword{ts}



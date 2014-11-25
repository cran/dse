\name{check.consistent.dimensions}
\alias{check.consistent.dimensions}
\alias{check.consistent.dimensions.default}
\alias{check.consistent.dimensions.ARMA}
\alias{check.consistent.dimensions.SS}
\alias{check.consistent.dimensions.TSdata}
\alias{check.consistent.dimensions.TSestModel}

\title{Check Consistent Dimensions}
\description{Check that dimensions of a model and data agree.}
\usage{
    check.consistent.dimensions(obj1, obj2=NULL)
    \method{check.consistent.dimensions}{default}(obj1, obj2=NULL)
    \method{check.consistent.dimensions}{ARMA}(obj1, obj2=NULL)
    \method{check.consistent.dimensions}{SS}(obj1, obj2=NULL)
    \method{check.consistent.dimensions}{TSdata}(obj1, obj2=NULL)
    \method{check.consistent.dimensions}{TSestModel}(obj1, obj2=NULL)
    }
\arguments{
    \item{obj1}{An object containing a TSmodel, TSdata, or TSestModel,
       depending on the method}
    \item{obj2}{Another object containing TSdata corresponding to the TSmodel in
       obj1, or a TSmodel corresponding to the TSdata in obj1.}
}
\value{logical}
\details{
    Check that dimensions of a model and data agree. If \code{obj1} is a
    \code{TSestModel} then if \code{obj2} is \code{NULL, TSdata} is 
    taken from \code{obj1}. 
}

\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    check.consistent.dimensions(model)
}
%\keyword{DSE}
\keyword{ts}


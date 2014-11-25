\name{coef.TSmodel.ee}
\alias{coef.TSmodel.ee}
\alias{coef.TSestModel.ee}
\title{Specific Methods for coef}
\description{See the generic function description.}
\usage{
    \method{coef}{TSmodel.ee}(object, criterion.args=NULL, ...)
    \method{coef}{TSestModel.ee}(object, criterion.args=NULL, ...)
}
\arguments{
    \item{object}{an object (model) from which to extract
       coefficients(parameters).}
    \item{criterion.args}{arguments to be passed to this method when it is
         called by \code{EstEval}.}
   \item{...}{(further arguments, currently disregarded).}
}
\details{
    The methods ***.ee are intended mainly to be called from \code{EstEval} as
    criterion for evaluating an estimation method. See \code{coef}. }
\seealso{
    \code{\link{EstEval}} 
    \code{\link{coef}} }
%\keyword{DSE}
\keyword{ts}


\name{roots.estimated.models}
\alias{roots.estimated.models}
\alias{roots.forecastCovEstimatorsWRTtrue}
\alias{roots.coef.ee}
\alias{roots.roots.ee}
\alias{roots.TSestModel.ee}
\alias{roots.TSmodel.ee}

\title{Roots Specific Methods}
\description{See the generic function description.}
\usage{
    \method{roots}{estimated.models}(obj, digits=options()$digits, mod =FALSE)
    \method{roots}{forecastCovEstimatorsWRTtrue}(obj, digits=options()$digits,
       mod=FALSE)
    \method{roots}{coef.ee}(obj, criterion.args=NULL)
    \method{roots}{roots.ee}(obj, ...)
    \method{roots}{TSestModel.ee}(obj, criterion.args=NULL)
    \method{roots}{TSmodel.ee}(obj, criterion.args=list(randomize = TRUE))
}
\arguments{
    \item{obj}{an object from which roots are to be extracted or calculated and
       printed.}
    \item{digits}{an integer indicating the number of significant digits to be
       printed (passed to the print method).}
    \item{mod}{if TRUE the modulus of the roots is calculated. Otherwise, a
       complex value may result.}
    \item{criterion.args}{arguments to be passed to this method when it is
         called by \code{EstEval}.}
    \item{...}{arguments to be passed to other methods.}
}
\details{
    The methods ***.ee are intended mainly to be called from EstEval as
    criterion for evaluating an estimation method.
}
\seealso{\code{\link{roots}}
         \code{\link{stability}}
         \code{\link{EstEval}}
	}

%\keyword{DSE}
\keyword{ts}


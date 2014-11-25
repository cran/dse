\name{roots.estimatedModels}
\alias{roots.estimatedModels}
\alias{roots.forecastCovEstimatorsWRTtrue}
\alias{roots.coefEstEval}
\alias{roots.rootsEstEval}
\alias{roots.TSestModelEstEval}
\alias{roots.TSmodelEstEval}

\title{Roots Specific Methods}
\description{See the generic function description.}
\usage{
    \method{roots}{estimatedModels}(obj, digits=options()$digits, mod =FALSE)
    \method{roots}{forecastCovEstimatorsWRTtrue}(obj, digits=options()$digits,
       mod=FALSE)
    \method{roots}{coefEstEval}(obj, criterion.args=NULL)
    \method{roots}{rootsEstEval}(obj, ...)
    \method{roots}{TSestModelEstEval}(obj, criterion.args=NULL)
    \method{roots}{TSmodelEstEval}(obj, criterion.args=list(randomize = TRUE))
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
\seealso{\code{\link[dse1]{roots}}
         \code{\link[dse1]{stability}}
         \code{\link{EstEval}}
	}

%\keyword{DSE}
\keyword{ts}


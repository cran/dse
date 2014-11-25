\name{tfplot.roots.ee}
\alias{tfplot.roots.ee}
\alias{plot.roots.ee}

\title{Specific tfplot methods for roots.ee (EstEval) objects}
\description{See the generic function description.}
\usage{
    \method{tfplot}{roots.ee}(x, ...)
    \method{plot}{roots.ee}(x, complex.plane=TRUE, cum=TRUE, norm=FALSE,
       bounds=TRUE, transform=NULL, invert=FALSE, Sort=TRUE, ...)
}
\arguments{
    \item{x}{an object for which a tfplot is to be produced.}
    \item{complex.plane}{logical indicating if the plot should be on the complex
       plane.}
    \item{cum}{logical indicating if the cumulative sum of roots should be
    plotted}
    \item{invert}{logical indicating if the inverse of roots should be
    plotted}
    \item{Sort}{logical indicating if the roots should be sorted.}
    \item{...}{arguments passed to other methods.}
    \item{norm}{logical indicating if the euclidean norm of roots should be
       plotted (square root of the sum of squared roots).}
    \item{bounds}{logical indicating if estimated one standard error bounds
       should be plotted around the lines for the true roots.}
    \item{transform}{an optional string indicating the name of a function which
       should be applied to the roots before plotting.}
}

\details{
If complex.plane is TRUE then all results are plotted on a complex plane and 
   the arguements cum and Sort do not apply. If complex.plane is FALSE 
   then a sequential plot of the real and imaginary parts is produced.
 If cum is true the cummulative average is plotted.
 If mod is true the modulus is used, otherwise real and imaginary are separated.
 if invert is true the reciprical is used (before cummulating).
 if Sort is true then sort is applied (before cum but after mod) by the Re part of the root.
   Some grouping is usually necessary since roots are not in an obvious order
   but sorting by the real part of the roots could be improved upon.
}
\seealso{
    \code{\link{tfplot}}
    \code{\link{EstEval}}
}

%\keyword{DSE}
\keyword{ts}


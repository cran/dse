\name{tfplot.coef.ee}
\alias{tfplot.coef.ee}

\title{Specific tfplot methods for coef.ee (EstEval) objects}
\description{See the generic function description.}
\usage{
    \method{tfplot}{coef.ee}(x, cum=TRUE, norm=FALSE, bounds=TRUE,
        invert=FALSE, Sort=FALSE, graphs.per.page = 5)
}
\arguments{
    \item{x}{an object for which a tfplot is to be produced.}
    \item{cum}{logical indicating if the cumulative sum of roots should be
       plotted}
    \item{invert}{logical indicating if the inverse of roots should be
       plotted}
    \item{Sort}{logical indicating if the roots should be sorted.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{norm}{logical indicating if the euclidean norm of roots should be
       plotted (square root of the sum of squared roots).}
    \item{bounds}{logical indicating if estimated one standard error bounds
       should be plotted around the lines for the true roots.}
}

\details{
 If cum is true the cummulative average is plotted.
 If norm is true the norm is used, each parameter is plotted.
 If invert is true the reciprical is used (before cummulating).
 If Sort is true then sort is applied (before cum). This is not usually
   recommended but of interest
   with estimation methods like black.box which may not return parameters
   of the same length or in the same order.
 Plotting the true lines only makes sense if truth is the same length as 
  result (and sometimes not even then). 
}
\seealso{
    \code{\link{tfplot}}
    \code{\link{EstEval}}
}
%\keyword{DSE}
\keyword{ts}


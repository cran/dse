\name{00Intro.setRNG}
\alias{00Intro.setRNG}
\title{setRNG}

\section{Usage}{
   library("setRNG") 
}

\description{
  Programs to set random number generator (and seed) in R and S.}
\section{Introduction}{
This library provides tools to reset the (normal) random number generator to a
reproducible state. It uses the R/S tools for setting the seed, but also sets
the mechanism for converting uniform numbers to normally distributed numbers.

The library also implements an approach to random number generation 
which allows the same random experiments to be replicated in S and R.
As of R 1.0.1 this adds nothing to R except for the function setRNG which sets
both RNGkind and the seed. The functions allow the R results using 
Wichmann-Hill and Box-Muller to be replicated in S.
These functions are intended primarily as a way to confirm that other code works in 
the same way in both S and R, not as an improved RNG. Default and other
RNGs can  still be used and are probably preferred for both speed and 
theoretical reasons. The present version of this works only for uniform 
and normal distributions. These were done with the aid of an example 
from B. D. Ripley.

The files in the S/ subdirectory of the package are for use with S (not R).
}
  

\examples{
is.R()
setRNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
rnorm(10)
}
\keyword{programming}
\keyword{interface}
\keyword{utilities}
\keyword{distribution}


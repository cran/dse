\name{Riccati}
\alias{Riccati}
\title{Riccati Equation}
\description{Solve a Matrix Riccati Equation}
\usage{
    Riccati(A, B, fuzz=1e-10, iterative=FALSE)
}
\arguments{
    \item{A}{A matrix.}
    \item{B}{A matrix.}
    \item{fuzz}{The tolerance used for testing convergence.}
    \item{iterative}{If TRUE an iterative solution technique is used.}
}
\value{xxx}
\details{
    Solve  Riccati equation  P = APA' + B 
    by eigenvalue decompostion of a symplectic matrix or by iteration.
}
\references{
    R. J. Vaccaro and T.Vukina(1993),
      "A Solution to the Positivity Problem in the State-Space Approach 
      to Modeling Vector-Valued Time Series." 
      \emph{Journal of Economic Dynamics and Control}, 17,
    Anderson & Moore, (1979) sec 6.7,  
    D. Vaughan (1970) \emph{IEEE Tr AC}.
    A. Laub (1983), \emph{Proc IEEE conf Decision and Control}. 
    Gudmundsson, Kenney and Laub (1992) \emph{IEEE Tr AC}.
}
\seealso{
    \code{\link{eigen}}
}
\note{This procedure has not been tested.}

\concept{DSE}
\keyword{ts}
\keyword{algebra}


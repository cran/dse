\name{Riccati}
\alias{Riccati}
\title{Riccati Equation}
\description{Solve a Matrix Riccati Equation}
\usage{
    Riccati(A, B, fuzz=1e-10, iterative=F)
}
\arguments {
    \item{A}{A matrix.}
    \item{B}{A matrix.}
    \item{fuzz}{The tolerance used for testing convergence.}
    \item{iterative}{If T an iterative solution technique is used.}
}
\value{xxx}
\details{
    Solve  Riccati equation  P = APA' + B 
    by eigenvalue decompostion of a symplectic matrix or by iteration.
}
\references{
    Vaccaro and Vukina, 
    Anderson & Moore, 1979 sec 6.7,  
    D.Vaughan (Oct 1970) \emph{IEEE Tr AC}.
    A. Laub (1983), \emph{Proc IEEE conf Decision and Control}. 
      Gudmundsson, Kenney and Laub (Apr 1992) \emph{IEEE Tr AC}.
}
\seealso{
    \code{\link{eigen}}
}
\note{This procedure has not been tested.}

%\keyword{DSE}
\keyword{ts}
\keyword{algebra}


\name{SS}
\alias{SS}
\alias{is.SS}
\alias{is.innov.SS}
\alias{is.non.innov.SS}

\title{State Space Models}
\description{Construct a }
\usage{
    SS(F.=NULL, G=NULL, H=NULL, K=NULL, Q=NULL, R=NULL, z0=NULL, P0=NULL,
             description=NULL, names=NULL, input.names=NULL, output.names=NULL)
    is.SS(obj)
    is.innov.SS(obj)
    is.non.innov.SS(obj)
}
\arguments{
    \item{F.}{(nxn) is the state transition matrix F.}
    \item{H}{(pxn)is the output matrix H.}
    \item{Q}{
    (nxn) is the input matrix of the system noise and the noise is 
    assumed to be white. Some authors (eg. Harvey) modify this 
    as rt*qt*rt' where rt is the matrix for the system noise and qt is 
    the noise cov, but that is redundant.
    }
    \item{R}{
    (pxp) is the input matrix of the output (measurement) noise, which 
    is assumed white. (probably need R if p>n )
    }
    \item{G}{(nxp)is the control (input) matrix.
        G should be NULL if there is no input.}
    \item{K}{(nxp)is the Kalman gain.}
    \item{z0}{vector indicating estimate of the state at time 0.
       Set to zero if not supplied.}
    \item{P0}{a matrix indicating initial tracking error P(t=1|t=0).
       Set to I if not supplied.}
    \item{description}{String. An arbitrary description.}
    \item{names}{A list with elements input and output, each a vector of 
       strings. Arguments input.names and output.names should not be used if
       argument names is used.}
    \item{input.names}{
      A vector of character strings indicating input variable names.
    }
    \item{output.names}{
      A vector of character strings indicating output variable names.
    }
    \item{obj}{an object.}
}
\value{An SS TSmodel}
\details{
State space models have a
further sub-class: innov or non-innov, indicating an innovations form
or a non-innovations form.  

The state space (SS) model is defined by:

      z(t) = Fz(t-1) + Gu(t) + Qe(t)
      y(t) = Hz(t)  + Rw(t)

or the innovations model:

       z(t) = Fz(t-1) + Gu(t) + Kw(t-1)
       y(t) = Hz(t)  + w(t)

    \item{F}{(nxn) is the state transition matrix F.}
    \item{H}{(pxn)is the output matrix H.}
    \item{Q}{
    (nxn) is the input matrix of the system noise and the noise is 
    assumed to be white. Some authors (eg. Harvey) modify this 
    as rt*qt*rt' where rt is the matrix for the system noise and qt is 
    the noise cov, but that is redundant.
}
    \item{R}{
    (pxp) is the input matrix of the output (measurement) noise, which 
    is assumed white. (probably need R if p>n )
    }
    \item{G}{(nxp)is the control (input) matrix.}
    \item{K}{(nxp)is the Kalman gain.}
    \item{y}{is the p dimensional output data.}
    \item{u}{is the m dimensional exogenous (input) data.}
    \item{z}{
    is the n dimensional (estimated) state at time t,  E[z(t)|y(t-1), u(t)]
    denoted E[z(t)|t-1]. An initial value for z can be specified as z0 and 
    an initial one step ahead state tracking error (for non-innovations 
    models) as P0.
}
    \item{z0}{An initial value for z can be specified as z0.}
    \item{P0}{
    An initial one step ahead state tracking error (for non-innovations 
    models) can be specified as P0.
}  
    \item{K, Q, R}{
    For sub-class innov the Kalman gain K is specified but not Q and R.
    For sub-class non-innov Q and R are specified but not the Kalman gain K.
}
}

\seealso{
    \code{\link{TSmodel}}
    \code{\link{ARMA}}
}
\examples{
    f <- array(c(.5,.3,.2,.4),c(2,2))
    h <- array(c(1,0,0,1),c(2,2))
    k <- array(c(.5,.3,.2,.4),c(2,2))
    ss <- SS(F=f,G=NULL,H=h,K=k)
    is.SS(ss)
}

%\keyword{DSE}
\keyword{ts}


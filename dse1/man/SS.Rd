\name{SS}
\alias{SS}
\alias{is.SS}
\alias{is.innov.SS}
\alias{is.nonInnov.SS}

\title{State Space Models}
\description{Construct a }
\usage{
    SS(F.=NULL, G=NULL, H=NULL, K=NULL, Q=NULL, R=NULL, z0=NULL, P0=NULL,
             description=NULL, names=NULL, input.names=NULL, output.names=NULL)
    is.SS(obj)
    is.innov.SS(obj)
    is.nonInnov.SS(obj)
}
\arguments{
    \item{F.}{(nxn) state transition matrix.}
    \item{H}{(pxn) output matrix.}
    \item{Q}{(nxn) matrix specifying the system noise distribution.}
    \item{R}{(pxp) matrix specifying the output (measurement) noise distribution.}
    \item{G}{(nxp) input (control) matrix. G should be NULL if there is no input.}
    \item{K}{(nxp) matrix specifying the Kalman gain.}
    \item{z0}{vector indicating estimate of the state at time 0.
       Set to zero if not supplied.}
    \item{P0}{matrix indicating initial tracking error P(t=1|t=0).
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

      z(t) =Fz(t-1) + Gu(t) + Qe(t)
      
      y(t) = Hz(t)  + Rw(t)

or the innovations model:

       z(t) =Fz(t-1) + Gu(t) + Kw(t-1)
       
       y(t) = Hz(t)  + w(t)

Matrices are as specified above in the aguments, and 
    \item{y}{is the p dimensional output data.}
    \item{u}{is the m dimensional exogenous (input) data.}
    \item{z}{is the n dimensional (estimated) state at time t,  
    E[z(t)|y(t-1), u(t)] denoted E[z(t)|t-1]. An initial value for z can 
    be specified as z0 and an initial one step ahead state tracking 
    error (for non-innovations models) as P0.}
    \item{z0}{An initial value for z can be specified as z0.}
    \item{P0}{An initial one step ahead state tracking error (for 
    non-innovations models) can be specified as P0.}  
    \item{K, Q, R}{
    For sub-class innov the Kalman gain K is specified but not Q and R.
    For sub-class non-innov Q and R are specified but not the Kalman gain K.}
    \item{e and w}{are typically assumed to be white noise in the 
    non-innovations form, in which case
    the covariance of the system noise is QQ' and the covariance of 
    the measurement noise is RR'. The covariance of e and w  can be specified 
    otherwise in the simulate 
    method \code{simulate.SS} for this class of model, but the assumption is
    usually maintained when estimating models of this form (although, not by all
    authors).
    
    Typically, an non-innovations form is harder to identify than an innovations
    form. Non-innovations form would typically be choosen when there is
    considerable theoretical or physical knowledge of the system (e.g. the
    system was built from known components with measured physical values).
    
    By default, elements in parameter matrices are treated as constants if they
    are exactly 1.0 or 0.0, and as parameters otherwise. A value of 1.001 would
    be treated as a parameter, and this is the easiest way to initialize an
    element which is not to be treated as a constant of value 1.0. Any matrix
    elements can be fixed to constants of any value using fixConstants.  
    }
}

\reference{B. D. O.Anderson & J. B.Moore \emph{Optimal Filtering} p.39,44.}

\seealso{
    \code{\link{TSmodel}}
    \code{\link{ARMA}}
    \code{\link{simulate.SS}}
    \code{\link{l.SS}}
    \code{\link{fixConstants}}
}
\examples{
    f <- array(c(.5,.3,.2,.4),c(2,2))
    h <- array(c(1,0,0,1),c(2,2))
    k <- array(c(.5,.3,.2,.4),c(2,2))
    ss <- SS(F=f,G=NULL,H=h,K=k)
    is.SS(ss)
    ss
}

%\keyword{DSE}
\keyword{ts}


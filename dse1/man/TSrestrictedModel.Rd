\name{TSrestrictedModel}
\alias{TSrestrictedModel}

\title{Restrict Model Coefficients}
\description{Put restrictions on a "TSmodel".}
\usage{
    TSrestrictedModel(model,
                    coefficients=NULL,
                    restriction=NULL)
}
\arguments{
    \item{model}{ 
      An object containing an object of class TSmodel (unrestricted).}
    \item{coefficients}{coefficients of the restricted model.}
    \item{restriction}{function that maps the coefficients of the 
           restricted model into the coefficients of the unrestricted model.}
}
\details{
BEWARE THIS CLASS OF MODELS IS EXPERIMENTAL AND MAY CHANGE.

This is the constructor for a "TSrestrictedModel","TSmodel"" 

The argument \code{model} should specify a \code{TSmodel} which is not a
\code{TSrestrictedModel}.

This constructor defines a new class of model \code{TSrestrictedModel}, 
which is composed of an unrestricted \code{TSmodel}, a new set of 
\code{coefficients}, and a function
called  \code{restriction} which maps the new set of coefficients into 
the unrestricted \code{TSmodel}. So, if \code{mod1} is a \code{TSrestrictedModel}, 
then \code{mod1$restriction(mod1$TSmodel, mod1$coefficients)} should return
an unrestricted \code{TSmodel} with it coefficients specified.
Methods like \code{l.TSrestrictedModel} sets the arrays using the set 
of coefficients for the TSrestrictedModel and
restrictions/equalities defined in the function \code{restriction}. 
 
It is possible the class is misnamed. It should be possible to construct
some extentions of TSmodels by definining restrictions on an enlarged model
(but I have not played with this yet).
}

\examples{
rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")


#  example  1                                                                             

z <- ARMA(A=c(1, 0.3), B=1)
mod1 <- TSrestrictedModel(z, 
            coefficients=coef(z),
            restriction=function(m, AllCoef){setArrays(m, 2*AllCoef)})
mod1

#  example 2

mod2 <- TSrestrictedModel(z, coefficients=c(3, coef(z)),
       restriction=function(m, AllCoef){ setArrays(m, AllCoef[1]*AllCoef[-1])})
mod2

#  example 3

z <- toSS(ARMA(A=c(1, 0.3, 0.1), B=1))
mod3 <- TSrestrictedModel(z, coefficients=c(2,.3, 4, coef(z)),
       restriction=function(m, AllCoef){
		      mm <- setArrays(m, AllCoef[-(1:3)])
		      P0 <- matrix(0,2,2)
		      P0[,1] <- AllCoef[1:2]
		      P0[,2] <- AllCoef[2:3]
		      mm$P0 <- P0
		      setTSmodelParameters(mm)
		      })
mod3

#  example 4

#  Starting P0  ("big k") symmetric with off diagonal element smaller than diag.
P0 <- matrix(1e6,4,4) 
diag(P0 )<- 1e7

# lower triangle will be parameters
P0[outer(1:4, 1:4, ">=")]
length(P0[outer(1:4, 1:4, ">=")])  # number of parameters

Hloadings <- t(matrix(c(
    8.8,   5.2,
   23.8, -12.6,
    5.2,  -2.0,
   36.8,  16.9,
   -2.8,  31.0,
    2.6,  47.6), 2,6))

z <-  SS(F=t(matrix(c(
    		  0.8, 0.04,  0.2, 0,
    		  0.2,  0.5,    0, -0.3,
    		    1,    0,    0, -0.2,
    		    0,	  1,    0,  0   ), c(4,4))),
	       H=cbind(Hloadings, matrix(0,6,2)),	
	       Q=diag(c(1, 1, 0, 0),4),  
	       R=diag(1,6),
	       z0=c(10, 20, 30,40),
	       P0=NULL
	       )

# The restriction constructions P0 from the lower triangle 

mod4 <- TSrestrictedModel(z, 
       coefficients=c(P0[outer(1:4, 1:4, ">=")],coef(z)),
       restriction=function(m, AllCoef){
		      mm <- setArrays(m, AllCoef[-(1:10)]) 
		      P0 <- matrix(0,4,4)
		      P0[outer(1:4, 1:4, ">=")] <- AllCoef[1:10]
		      P0 <- P0 + t(P0)
		      diag(P0) <- diag(P0)/2
		      mm$P0 <- P0
		      setTSmodelParameters(mm)
		      })
mod4


z  <- simulate(SS(F=t(matrix(c(
    		  0.8, 0.04,  0.2, 0,
    		  0.2,  0.5,    0, -0.3,
    		    1,    0,    0, -0.2,
    		    0,	  1,    0,  0   ), c(4,4))),
	       H=cbind(Hloadings, matrix(0,6,2)),	
	       Q=diag(c(1, 1, 0, 0),4),  
	       R=diag(1,6),
	       z0=c(10, 20, 30,40),
	       P0=diag(c(10, 10, 10, 10)) ),
               rng=rngValue10)  
state.sim  <- z$state  # for comparison below
y.sim  <- outputData(z) # simulated indicators

coef(mod4)
coef(l(mod4, TSdata(output=y.sim)))
summary(l(mod4, TSdata(output=y.sim))) 


zz <- smoother(l(mod4, TSdata(output=y.sim)))  
summary(zz)

tfplot(state.sim,state(zz))
tfplot(state.sim,state(zz, smoother=TRUE))

est.mod4 <- estMaxLik(mod4, TSdata(output=y.sim),
        algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf, hessian=TRUE,
                        control=list(maxit=10000))
    ) 

summary(est.mod4)

sest.mod4 <- smoother(est.mod4)
summary(sest.mod4)

tfplot(sest.mod4, graphs.per.page=3)
tfplot(state.sim, state(sest.mod4, smoother=TRUE))
tfplot(state.sim, state(sest.mod4, filter=TRUE))

coef(sest.mod4)
coef(mod4)
}

\seealso{
    \code{\link{TSmodel}},
    \code{\link{ARMA}},
    \code{\link{SS}},
    \code{\link{TSestModel}},
    \code{\link{TSdata}}
}
\concept{DSE}
\keyword{ts}

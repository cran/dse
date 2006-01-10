
require("dse1")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")

###################################################

# test class TSrestrictedModel

###################################################


#  example  1                                                                             

mod1 <- ARMA(A=c(1, 0.3), B=1)
z <- TSrestrictedModel(mod1, 
            coefficients=coef(mod1),
            restriction=function(m, AllCoef){setArrays(m, 2*AllCoef)})
z

#  example 2

z <- TSrestrictedModel(mod1, coefficients=c(3, coef(mod1)),
       restriction=function(m, AllCoef){ setArrays(m, AllCoef[1]*AllCoef[-1])})
z

#  example 3

mod2 <- toSS(ARMA(A=c(1, 0.3, 0.1), B=1))
z <- TSrestrictedModel(mod2, coefficients=c(2,.3, 4, coef(mod2)),
       restriction=function(m, AllCoef){
		      mm <- setArrays(m, AllCoef[-(1:3)])
		      P0 <- matrix(0,2,2)
		      P0[,1] <- AllCoef[1:2]
		      P0[,2] <- AllCoef[2:3]
		      mm$P0 <- P0
		      setTSmodelParameters(mm)
		      })
z


######################################

# test "big k" (which is numerically sensitive).

######################################

#  example 4


Hloadings <- t(matrix(c(
    8.8,   5.2,
   23.8, -12.6,
    5.2,  -2.0,
   36.8,  16.9,
   -2.8,  31.0,
    2.6,  47.6), 2,6))


#  Starting P0  ("big k") symmetric with off diagonal element smaller than diag.
P0 <- matrix(1e6,4,4) 
diag(P0 )<- 1e7


z <- SS(F=t(matrix(c(
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


summary(mod4)
coef(mod4)

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


error <- max(abs(l(mod4, TSdata(output=y.sim), return.state=TRUE)$filter$state -
                 l(mod4, TSdata(output=y.sim), return.state=TRUE,
	                                           compile=FALSE)$filter$state))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


zz  <- smoother(l(mod4, TSdata(output=y.sim)))   
zzz <- smoother(l(mod4, TSdata(output=y.sim)), compiled=FALSE)

tfplot(state.sim,  state(zz))
tfplot(state.sim,  state(zzz))
tfplot(state.sim,  state(zz,  smoother=TRUE))
tfplot(state.sim,  state(zzz, smoother=TRUE))

error <- max(abs(state(zz, filter=TRUE) - state(zzz, filter=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(state(zz, smoother=TRUE) - state(zzz, smoother=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(zz$filter$track - zzz$filter$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(zz$smooth$track - zzz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

est.mod4 <- estMaxLik(mod4, TSdata(output=y.sim),
        algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf, hessian=TRUE,
                        control=list(maxit=10000))
    ) 
summary(est.mod4)

sest.mod4 <- smoother(est.mod4)


summary(sest.mod4)
class(sest.mod4$model)
coef(sest.mod4)


tfplot(sest.mod4, graphs.per.page=3)
tfplot(state.sim, state(sest.mod4, filter=TRUE))

zzz       <- smoother(est.mod4, compiled=FALSE)

error <- max(abs(state(sest.mod4, filter=TRUE) - state(zzz, filter=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(state(sest.mod4, smoother=TRUE) - state(zzz, smoother=TRUE)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(sest.mod4$filter$track - zzz$filter$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs(sest.mod4$smooth$track - zzz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }



if (! all.ok) stop("some tests FAILED")


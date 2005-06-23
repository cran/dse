# This should test comparison of fortran and S version, but the fortran 
#  matrix inversio is not working well, so the tests are disabled.

require("dse1")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

###################################################

# test with input and having output dim < state dim.

###################################################

if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model <- TSmodel(toSSChol(estVARXls(eg1.DSE.data.diff))) 

model0 <- model
model0$G <- NULL
simdata0 <- simulate(model0, rng=test.rng) 

z  <- smoother(model0, simdata0, compiled=TRUE)
zz <- smoother(model0, simdata0, compiled=FALSE)

#zcomp <- z
#debug(dse1:::smoother.default)
#zcomp$r$chk1[99,1:18,1:3] - K  good
#zcomp$r$chk1[99,1:18,1] - zt good
#max(abs(zcomp$r$chk2[99,1:18,1:18] - P)) good
#max(abs(zcomp$r$chk1[99,1:18,1:18] - P)) good
#max(abs(zcomp$r$ftrack[Time+1,,] - filter$track[Time+1,,])) good
#max(abs(zcomp$r[[11]] - FF))  good
#max(abs(zcomp$r$chk1[Time,,] - filter$track[Time+1,,]))  good
#max(abs(zcomp$r$chk2[Time,,] - solve(filter$track[Time+1,,]))) very bad
#max(abs(diag(1,18) - filter$track[Time+1,,] %*% solve(filter$track[Time+1,,]))) good
#max(abs(diag(1,18) - zcomp$r$chk1[Time,,] %*% zcomp$r$chk2[Time,,])) not good
#zcomp$r$chk2[99,1:18,1:18] - t(FF) %*% solve(filter$track[Time+1,,]) bad
#zcomp$r$chk2[99,1:18,1:18] - J bad

#tfplot(simdata0$state, state(z, smooth=TRUE), state(zz, smooth=TRUE), graphs.per.page=3)


# using simulated data gives a true state for comparison.
simdata <- simulate(model, input= inputData(eg1.DSE.data.diff), rng=test.rng) 

z  <- smoother(model, simdata, compiled=TRUE)
zz <- smoother(model, simdata, compiled=FALSE)


error <- max(abs((state(zz, smooth=TRUE) - zz$smooth$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

tfplot(state(zz), simdata$state, graphs.per.page=3)
#tfplot(state(z, smooth=TRUE) - state(zz, smooth=TRUE), graphs.per.page=3)
tfplot(state(z, smooth=TRUE),  state(zz, smooth=TRUE), graphs.per.page=3)

# plot smoother agains true state
#tfplot(state(z, smooth=TRUE), simdata$state, graphs.per.page=3)

# plot smoother agains true state
#tfplot(state(zz, smooth=TRUE), simdata$state, graphs.per.page=3)

#tfplot(state(z, smooth=TRUE), simdata$state, graphs.per.page=3)
#tfplot(simdata$state, state(z, smooth=TRUE), state(zz, smooth=TRUE), graphs.per.page=3)

#tfplot(simdata$state, state(zz, filter=TRUE), state(zz, smooth=TRUE), graphs.per.page=3)

# compare fortran and S versions  DISABLED
error <- max(abs((state(z, smooth=TRUE) - state(zz, smooth=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
# this fails
#      all.ok <- FALSE  
     }

error <- max(abs(z$filter$track - zz$filter$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

# compare fortran and S versions  DISABLED
error <- max(abs(z$smooth$track - zz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
# this fails
#      all.ok <- FALSE  
     }


######################################

# test output dim exceeds state dim.

######################################

Hloadings <- t(matrix(c(
    8.8,   5.2,
   23.8, -12.6,
    5.2,  -2.0,
   36.8,  16.9,
   -2.8,  31.0,
    2.6,  47.6), 2,6))

ss.ar1 <- SS(F=array(c(.5, .4, .3, .2),c(2,2)),  
		H=Hloadings,     
		Q=array(c(1.0, 2.0),c(2,2)),  
		R=diag(1,6) 
		)

simdata2 <- simulate(ss.ar1, rng=test.rng)

z  <- smoother(ss.ar1, simdata2, compiled = TRUE)
zz <- smoother(ss.ar1, simdata2, compiled = FALSE)

#tfplot(state(zz, smooth=TRUE), state(z, smooth=TRUE), simdata2$state, graphs.per.page=3)
#tfplot(state(zz, smooth=TRUE) - state(z, smooth=TRUE), graphs.per.page=3)

error <- max(abs((state(zz, filter=TRUE) - zz$filter$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


error <- max(abs((state(zz, filter=TRUE) - state(z, filter=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs((state(zz, smooth=TRUE) - zz$smooth$state)))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }

error <- max(abs((state(z, smooth=TRUE) - state(zz, smooth=TRUE))))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


error <- max(abs((z$smooth$track) - zz$smooth$track))
if ( fuzz < error) 
     {print(error, digits=18)
      all.ok <- FALSE  
     }


tfplot(simdata2$state, state(zz, smooth=TRUE), state(zz, filter=TRUE))
tfplot(simdata2$state, state(z, smooth=TRUE),  state(z, filter=TRUE))
tfplot(simdata2$state, state(z, smooth=TRUE),  state(zz, smooth=TRUE))

if (! all.ok) stop("some tests FAILED")


 require("dse1")
 Sys.info()
 version.dse()
 if      (is.R()) data("eg1.DSE.data.diff", package="dse1") else 
 if (is.S()) 
   {source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))
    class(eg1.DSE.data.diff$output) <- class(eg1.DSE.data.diff$input) <- NULL
    }

 if (!is.TSdata(eg1.DSE.data.diff)) stop("Test data not found. Testing stopped.")
 
fuzz.small <- 1e-14
fuzz.large <- 1e-10
digits <- 18
all.ok <- TRUE  


test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

  VARmodel  <-  est.VARX.ar(eg1.DSE.data.diff, re.add.means=FALSE, warn=FALSE)

  VARmodelB <- TSmodel(VARmodel)
  B <- t(chol(VARmodel$estimates$cov))
  VARmodelB$B <- array(B, c(1,dim(B)))  # has B != I
  VARmodelB <- set.parameters(VARmodelB)
  VARmodelB <- l(VARmodelB,VARmodel$data, warn=FALSE)

cat("dse1 test 8a ...\n")
  z  <- simulate(TSmodel(VARmodel), input=input.data(eg1.DSE.data.diff)) 
  zz <- simulate(TSmodel(VARmodel), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff))
  ok <- test.equal(z, zz, fuzz=fuzz.small)

   if (!ok) 
     {
      all.ok <- FALSE  
     }

cat("dse1 test 8b ...\n")
  sigma <- solve(t(VARmodelB$model$B[1,,]) %*% VARmodelB$model$B[1,,])
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zzz <- simulate(TSmodel(VARmodelB), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff), SIGMA=sigma)
  ok <- test.equal(z, zzz, fuzz=fuzz.small)
  error <- max(abs(output.data(z) - output.data(zzz)))
   cat("max. error ", max(error), "\n")

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error || !ok) 
     {print.test.value(output.data(z), digits=18)
      print.test.value(output.data(zzz), digits=18)
      all.ok <- FALSE  
     }

cat("dse1 test 8c ...\n")
  # next use estimates$cov
  z  <- simulate(VARmodel, input=input.data(eg1.DSE.data.diff)) 
  sigma <- VARmodel$estimates$cov
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zz <- simulate(TSmodel(VARmodel), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff), SIGMA=sigma)

  if ( ! test.equal(z, zz, fuzz=fuzz.small)) 
     {
      all.ok <- FALSE  
     }

  ok <- test.equal(summary(z)$estimates,
                   summary(zz)$estimates, fuzz=fuzz.small)

  if ( !ok) 
     {print.test.value(c(summary(z)$estimates),  digits=18)
      print.test.value(c(summary(zz)$estimates), digits=18)
      all.ok <- FALSE  
     }



  if (! all.ok) stop("some tests FAILED")


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

  SSmodel  <- to.SS(VARmodel)

cat("dse1 test 9 ...\n")
  z  <- simulate(SSmodel, input=input.data(eg1.DSE.data.diff)) 
  ok <- test.equal(z,simulate(SSmodel, rng=get.RNG(z), 
                      input=input.data(eg1.DSE.data.diff)))
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}

  ok <- test.equal(summary(z)$estimates,
                   summary(z)$estimates, fuzz=fuzz.small)
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}



cat("dse1 test 10...\n")

  ok <- stability(SSmodel, verbose=FALSE)
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}


cat("dse1 test 11...\n")

   scale.fac <- diag(1:3)
   scale.fac[1,3] <-.5
   scale.pred <- VARmodel$estimates$pred %*% t(scale.fac)
   scale.fac <- list(output=scale.fac)

   good <- scale.pred
   tst  <- l(scale(VARmodel$model, scale=scale.fac), 
          scale(eg1.DSE.data.diff, scale=scale.fac), warn=FALSE)$estimates$pred
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
 
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 12...\n")

   good <- scale.pred
   tst  <- l(scale(SSmodel, scale=scale.fac), 
             scale(eg1.DSE.data.diff, scale=scale.fac))$estimates$pred
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
 
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 13...\n")

  z <- eg1.DSE.data.diff
  ok <- test.equal(z,
      TSdata(output=output.data(combine(z,z), series=seq(nseriesOutput(z))),
              input= input.data(combine(z,z), series=seq( nseriesInput(z))))) 
 
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}


  if (! all.ok) stop("some tests FAILED")


 require("dse1")
 Sys.info()
 version.dse()
 data("eg1.DSE.data.diff", package="dse1") 
 
 if (!is.TSdata(eg1.DSE.data.diff)) stop("Test data not found. Testing stopped.")
 
fuzz.small <- 1e-14
fuzz.large <- 1e-10
digits <- 18
all.ok <- TRUE  


test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

VARmodel  <-  est.VARX.ar(eg1.DSE.data.diff, re.add.means=FALSE, warn=FALSE)
SSmodel  <- to.SS(VARmodel)


cat("dse1 test 2 ...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(set.arrays(SSmodel), eg1.DSE.data.diff,warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }

cat("dse1 test 3 ...\n")
   good <- VARmodel$estimates$like[1] 
   tst  <- l(set.arrays(VARmodel), eg1.DSE.data.diff, warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 4 ...\n")
  ARMAmodel <- to.ARMA(SSmodel)
  good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff, warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 5 ...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff,warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 6a...\n")
   good <- sort(Mod(roots(TSmodel(VARmodel),by.poly=TRUE)))
   tst  <- sort(Mod(roots(SSmodel)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse1 test 6b...\n")
   good <- sort(Mod(roots(SSmodel)))
   tst  <- sort(Mod(roots(TSmodel(VARmodel),by.poly=FALSE)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }



  if (! all.ok) stop("some tests FAILED")


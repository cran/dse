 require("dse1")
 Sys.info()
 version.dse()
 if      (is.R()) data("eg1.DSE.data.diff", package="dse1") else 
 if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

 if (!is.TSdata(eg1.DSE.data.diff)) stop("Test data not found. Testing stopped.")
 
fuzz.small <- 1e-14
fuzz.large <- 1e-10
digits <- 18
all.ok <- T  


test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

VARmodel  <-  est.VARX.ar(eg1.DSE.data.diff, re.add.means=F, warn=F)
SSmodel  <- to.SS(VARmodel)


cat("dse1 test 2 ...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(set.arrays(SSmodel), eg1.DSE.data.diff,warn=F)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("dse1 test 3 ...\n")
   good <- VARmodel$estimates$like[1] 
   tst  <- l(set.arrays(VARmodel), eg1.DSE.data.diff, warn=F)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }


cat("dse1 test 4 ...\n")
  ARMAmodel <- to.ARMA(SSmodel)
  good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff, warn=F)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }


cat("dse1 test 5 ...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff,warn=F)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }


cat("dse1 test 6a...\n")
   good <- sort(Mod(roots(TSmodel(VARmodel),by.poly=T)))
   tst  <- sort(Mod(roots(SSmodel)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }


cat("dse1 test 6b...\n")
   good <- sort(Mod(roots(SSmodel)))
   tst  <- sort(Mod(roots(TSmodel(VARmodel),by.poly=F)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }



  if (! all.ok) stop("some tests FAILED")


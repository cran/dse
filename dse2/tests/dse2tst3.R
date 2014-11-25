  require("mva"); require("ts"); require("dse2") # adds dse, tframe, and syskern
 #x11()
  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
             # width=6, height=8, pointsize=10,
   Sys.info()
   version.dse()
   random.number.test() 





dse3.function.tests <- function(verbose=T, synopsis=T, fuzz.small=1e-14, fuzz.large=1e-8, graphics=T)
{ max.error <- NA
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))


# The seed is not important for most of these tests, but AIC eliminates all
#  parameters occassionally in some model selection tests.
 test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
# set.RNG(test.rng)

  if (synopsis & !verbose) cat("All dse3 tests ...")
  if (verbose) cat("dse3 test 0 ... ")
  data <- eg1.DSE.data.diff
  input.data(data) <- NULL
  mod1 <- TSmodel(est.VARX.ls(data))
  mod2 <- TSmodel(est.VARX.ar(data, re.add.means=F, warn=F))
  ok <- is.TSmodel(mod1)
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 1 ... ")
  z <- monte.carlo.simulations(mod1, replications=5, quiet=T)
  ok <- is.monte.carlo.simulation(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 2 ... ")
  ok <- test.equal(z, monte.carlo.simulations(mod1, replications=5,
                                     rng=get.RNG(z), quiet=T))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 3 ... ")
  z <- eval.estimation(mod1, replications=3,  estimation="est.VARX.ls",
            estimation.args=NULL, criterion="TSmodel", quiet=T)
  ok <- is.estimation.evaluation(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 4 ... ")
  zz <-summary(parms(z), verbose=F)
  ok <- T   # could be improved
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 5 ... ")
  zz <- summary(roots(z), verbose=F)
  ok <- T   # could be improved
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 6 ... ")
  zz <- parms(z)
  ok <- is.estimation.evaluation(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 7 ... ")
  zz <- roots(z)
  ok <- is.estimation.evaluation(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 8a... ")
  z <- horizon.forecasts(mod1, data, horizons=c(6,12), discard.before=20)
  error <- max(abs( c(z$horizon.forecasts[,100,])  -
 c(0.0048425425521641824594, 0.0031489473295282835973, 0.0037730234730729999594,
 0.0024354234760485438289, 0.0040593859721713481878, 0.0031982930612152113414)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 8b... ")
  z <- horizon.forecasts(l(to.SS(mod1), data),
                         horizons=c(6,12), discard.before=20)
  error <- max(abs( c(z$horizon.forecasts[,100,]) -
 c(0.0048425425521641824594, 0.0031489473295282844646, 0.0037730234730729995257,
 0.0024354234760485446963, 0.0040593859721713499225, 0.0031982930612152122088)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse3 test 9 ... ")
  zzz<-l(mod1,simulate(mod1))
  zz<-forecast.cov(zzz, discard.before=50, horizons=1:4)
  ok <- is.forecast.cov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 10... ")
  ok <- test.equal(zz, forecast.cov(zzz$model, 
             data=zzz$data, discard.before=50, horizons=1:4))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 11... ")
  zz <-forecast.cov(mod1,mod2, data=data, discard.before=30, zero=T, trend=T)

  ok <- is.forecast.cov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 12... ")
  zzz <-forecast.cov(to.SS(mod1),to.SS(mod2), data=data, 
                 discard.before=30, zero=T, trend=T)

  ok <- test.equal(zz,zzz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 13... ")
  zz <-out.of.sample.forecast.cov.estimators.wrt.data(data,
               estimation.methods = list(est.VARX.ar= list(max.lag=2, warn=F), 
                                         est.VARX.ls= list(max.lag=2)))
  ok <- is.forecast.cov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 14... ")
  zz <- forecast.cov.wrt.true(list(mod1,mod2),mod1, 
               pred.replications=2, Spawn=F, quiet=T, trend=NULL, zero=T)
  ok <- is.forecast.cov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 15... ")
  ok <- test.equal(zz, forecast.cov.wrt.true(list(mod1,mod2),mod1, 
          pred.replications=2, Spawn=.SPAWN, quiet=T, trend=NULL, zero=T,
          rng=get.RNG(zz)))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 16... ")
#  zz <- forecast.cov.estimators.wrt.true(mod1, Spawn=.SPAWN, quiet=T, 
# This seems to cause a problem in Splus when .SPAWN is T although it may
#  work with default rng (at least it used to) but test.rng is now set
#  to give same results as in R.
  zz <- forecast.cov.estimators.wrt.true(mod1, Spawn=F, quiet=T, 
         estimation.methods=list(est.VARX.ls=NULL, est.VARX.ar=list(warn=F)), 
         est.replications=2, pred.replications=2, rng=test.rng)
  ok <- is.forecast.cov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 17... ")
  ok <- test.equal(zz, forecast.cov.estimators.wrt.true(mod1, Spawn=F, 
           estimation.methods=list(est.VARX.ls=NULL,est.VARX.ar=list(warn=F)), 
           est.replications=2, pred.replications=2, rng=get.RNG(zz)))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (graphics)
      {ok <- dse3.graphics.tests(verbose=verbose,  pause=F)
       all.ok <- all.ok & ok 
       if (verbose) cat("dse3 test 18 (graphics) ... ")
       if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
      }

  if (synopsis) 
    {if (verbose) cat("All dse3 tests completed")
     if (all.ok) cat(" OK\n") 
     else  cat(" some FAILED! max.error = ", max.error,"\n")
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}




dse3.graphics.tests <- function(verbose=T, synopsis=T)
{ if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (synopsis & !verbose) cat("dse3 graphics tests ...")
  if (verbose) cat("  dse3 graphics test 1 ...")
  # If no device is active then write to postscript file 
  if ( dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",
                   width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
            {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

# The seed is not important for most of these tests, but AIC eliminates all
#  parameters occassionally in some model selection tests.
 test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))

  data <- eg1.DSE.data.diff
  input.data(data) <- NULL
  output.data(data) <- output.data(data, series=1)  # [,1,drop=F]
  mod1 <- TSmodel(est.VARX.ls(data,max.lag=3))
  mod2 <- TSmodel(est.VARX.ar(data,max.lag=3, aic=F, warn=F))

  z <- eval.estimation(mod1, replications=10,  estimation="est.VARX.ls",
            estimation.args=list(max.lag=3), criterion="TSmodel", quiet=T)
  distribution(parms(z)) 
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 2 ...")
  distribution(roots(z))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 3 ...")
  z <- horizon.forecasts(mod1, data, horizons=c(6,12), discard.before=20)
  tfplot(z, start.=c(1985,1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 4 ...")
  zz <-forecast.cov(mod1,mod2, data=data,
                     discard.before=10, zero=T, trend=T)
  tfplot(zz)
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 5 ...")
  tfplot(zz, select.cov=c(1), select.trend=F)
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 6 ...")

  data <- eg1.DSE.data.diff
  input.data(data) <- NULL
# next causes Error ... all lags eliminated by AIC order selection.
#  output.data(data) <- output.data(data, series=1)  # [,1,drop=F]
  zz <-out.of.sample.forecast.cov.estimators.wrt.data(data,
               estimation.methods = list(est.VARX.ar=list(max.lag=2,warn=F),
                                         est.VARX.ls=list(max.lag=2))) 
  tfplot(zz, select.series=c(1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 7 ...")
  zz <- forecast.cov.wrt.true(list(mod1,mod2),mod1, rng=test.rng,
               pred.replications=2, Spawn=.SPAWN, trend=NULL, zero=T, quiet=T)
  tfplot(zz, select.cov=c(1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 8 ...")
  zz <- forecast.cov.estimators.wrt.true(mod1, Spawn=.SPAWN,  rng=test.rng,
         estimation.methods=list(est.VARX.ls=NULL,est.VARX.ls=list(max.lag=2)), 
#        estimation.methods=list(est.VARX.ls=NULL,est.VARX.ar=NULL), 
         est.replications=2, pred.replications=2, quiet=T)
  tfplot(zz, select.cov=c(1))
  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All dse3 graphics tests completed\n")
     else cat("completed\n")
    }
      
  invisible(T)
}



   dse3.function.tests(verbose=T, graphics=F) 
   dse3.graphics.tests(verbose=T)

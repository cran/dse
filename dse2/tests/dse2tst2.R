  require("mva"); require("ts"); require("dse2") # adds dse, tframe, and syskern
 #x11()
  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
             # width=6, height=8, pointsize=10,
   Sys.info()
   version.dse()
   random.number.test() 




dse2.function.tests <- function(verbose=T, synopsis=T, fuzz.small=1e-14, fuzz.large=1e-8, graphics=T)
{max.error <- NA
 if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
 else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

 if (synopsis & !verbose) cat("All dse2 tests ...") 
 if (verbose) cat("dse2 test 0 ... ")
  z <- eg1.DSE.data.diff
  z$input <- NULL
  mod1 <- TSmodel(est.VARX.ar(z, re.add.means=F, warn=F))
  ok <- is.TSmodel(mod1)
  all.ok <- ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse2 test 1 ... ")
  z <- est.black.box1(eg1.DSE.data.diff, verbose=F, max.lag=2)
  error <- max(abs(z$estimates$like[1]+4025.943051342767))
  ok <- is.TSestModel(z) &  (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("dse2 test 2 ... ")
  z <- est.wt.variables(eg1.DSE.data.diff, c(1,10,10),
                        estimation="est.VARX.ls")
  error <- max(abs(z$estimates$like[1]+4125.05572604540066)) 
  ok <- is.TSestModel(z) &   (fuzz.large > error)
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("dse2 test 3 ... ")
  z <- est.SS.Mittnik(eg1.DSE.data.diff, max.lag=2, n=3)
  error <- max(abs(z$estimates$like[1]+3794.0394069904219))
  ok <- is.SS(z$model) &   (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("dse2 test 4 ... ")
  z <- l( reduction.Mittnik(z, criterion="taic", verbose=F), 
         eg1.DSE.data.diff)
  error <- max(abs(z$estimates$like[1]+3795.6760513068380)) 
  ok <- is.SS(z$model)  &  (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  modSS<-z

  if (verbose) cat("dse2 test 5 ... ")
  z <- featherForecasts( modSS,  from.periods=c(250,300))
  error <- max(abs
       (c(z$featherForecasts[[1]][286,],z$featherForecasts[[2]][336,])
       -c(-0.00092229286770808757701, -0.0086020067525247358164, 
           0.0043454851777852505565,  -0.0066741302949233430319,
          -0.0089398331205012854933,   0.0021769124280658046222)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("dse2 test 6 ... ")
  # previously end=c(1969,6) when .diff data had wrong start date
  output.data(modSS$data) <- tfwindow(output.data(modSS), end=c(1969,7))
  # it should be possible to do the following instead, but tsp seems to
  # sometimes get mixed up in forecast and cause System terminating: bad address
  # output.data(modSS$data) <- output.data(modSS$data)[1:100,]
  z <- forecast(modSS, percent=c(90,100,110))

# previously 136 below
  error <- max(abs(
    c(z$forecast[[1]][36,],z$forecast[[2]][36,], z$forecast[[3]][36,])
     -c(-0.00310702417651131587, -0.00604105559321206804,0.00214657444656118738,
      -0.00345224972784219028, -0.00671228396225603124,0.00238508249578931863,
      -0.00379747527917305948, -0.00738351233129999531,0.00262359054501745074)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

if (graphics) 
 {if (verbose) cat("dse2 test 7 (graphics) ... ")
  ok <- dse2.graphics.tests(verbose=verbose, pause=T)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
 }

  if (synopsis) 
    {if (verbose) cat("All dse2 tests completed")
     if (all.ok) cat(" OK\n")
     else    
       {cat(", some FAILED!")
        if(max.error > fuzz.small)
            cat(" max. error magnitude= ", max.error,")")
        cat("\n")
       }
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}

dse2.graphics.tests <- function(verbose=T, synopsis=T)
{# graphics tests do not do any value comparisons
  if (synopsis & !verbose) cat("dse2 graphics tests ...")
  
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (verbose) cat("  dse2 graphics test 1 ...")

  # If no device is active then write to postscript file 
  if (dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

  data <- eg1.DSE.data.diff
  mod1 <- TSmodel(est.VARX.ls(data,max.lag=3))
  modSS <- l(to.SS(mod1),data)

  z <- featherForecasts( modSS,  from.periods=c(230,250))
  tfplot(z, start.=c(1980,1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse2 graphics test 2 ...")
  z <- forecast(modSS, percent=c(90,100,110))
  tfplot(z, start.=c(1985,1))
  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All dse2 graphics tests completed\n")
     else cat("completed\n")
    }
      

  invisible(T)
}



   dse2.function.tests(verbose=T, graphics=F)  
   dse2.graphics.tests(verbose=T)

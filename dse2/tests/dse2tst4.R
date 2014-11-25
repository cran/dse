  require("mva"); require("ts"); require("dse2") # adds dse, tframe, and syskern
 #x11()
  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
             # width=6, height=8, pointsize=10,
   Sys.info()
   version.dse()
   random.number.test() 




dse4.function.tests <- function(verbose=T, synopsis=T, 
		fuzz.small=1e-14, fuzz.large=1e-7, graphics=T)
{max.error <- 0
 if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
 if (is.S()) 
   {source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))
    class(eg1.DSE.data.diff$output) <- class(eg1.DSE.data.diff$input) <- NULL
    }

 if (synopsis & !verbose) cat("All dse4 tests ...") 
 if (verbose) cat("dse4 test 1 ... ")
  z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
  ok <- is.forecastCov.estimatorsWRTdata.subsets(z)
  all.ok <-  ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse4 test 2 ... ")
  z1 <- z$multi.model[[
       select.forecastCov(z, select.cov.best=1, verbose=F)$selection.index[2]]]
  subdata <- TSdata(output=output.data(eg1.DSE.data.diff, series=1:3))
  z2 <- estimate.models(subdata, estimation.sample =182, quiet = T, 
           estimation.methods = list(est.VARX.ls=list(max.lag=3)))
  output.data(subdata) <- output.data(subdata)[1:182,,drop=F]
#  input.data(subdata)  <- input.data(subdata) [1:182,,drop=F] not in subdata
  z3 <- est.VARX.ls(subdata, max.lag=3)
  ok <-      test.equal(z2$multi.model[[1]],z3$model)
  ok <- ok & test.equal(z2$multi.model[[1]],  z1)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse4 test 3 ... ")
#Rbug needs stepwise
if (is.R()) warning("skipping test 4 (requires stepwise).")
else
 {
   all.data <- TSdata(input=eg1.DSE.data.diff$output, 
                   output=eg1.DSE.data.diff$input )
   umodel <- build.input.models(all.data, max.lag=2)
   umodel <- build.diagonal.model(umodel)
   z  <- TSdata(output=output.data(all.data), 
                input=input.data(all.data, series=1:2))
	# previously ??input=input.data(extract(all.data, outputs=1, inputs=1:2)))
  ymodel <- est.VARX.ls(z, max.lag=3)$model 
  z <- ymodel$C
  ymodel$C <- array(0, c(dim(z)[1:2], nseriesOutput(umodel))) 
  ymodel$C[1:(dim(z)[1]), 1:(dim(z)[2]), 1:(dim(z)[3])] <- z 
  sim.data <- gen.mine.data(umodel, ymodel,
    rng= list(kind="default",seed=c(21,46,16,12, 51, 2, 31, 8, 42, 60, 7, 3)) )
  m.step <- mine.stepwise(sim.data, method="backward", plot.=F)
  error <- max(abs(m.step$stepwise$rss[c(1,27)] -
                c(5.65168201726030333e+01, 3.06441399376576440e+06)))
  # previously  c(47.537312899054931847,   4088283.2706551752053)))
  ok <- fuzz.large > error
  if (!ok) max.error <- max(error, max.error)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
 }

  if (graphics)
   {ok <- dse4.graphics.tests(verbose=verbose, pause=F)
    all.ok <- all.ok & ok 
   }

  if (synopsis) 
    {if (verbose) cat("All dse4 tests completed")
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

dse4.graphics.tests <- function(verbose=T, synopsis=T)
{ if (synopsis & !verbose) cat("dse4 graphics tests ...")
  if (verbose) cat("  dse4 graphics test 1 ...")

  # If no device is active then write to postscript file 
  if ( dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

  z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
  zz <- tfplot(z)

  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All dse4 graphics tests completed\n")
     else cat("completed\n")
    }
      
  invisible(T)
}



   dse4.function.tests(verbose=T, graphics=F) 
   dse4.graphics.tests(verbose=T)  #     test 3 needs stepwise

#   2000/04/25 10:03:23
# For installation instructions see the file read.me or the brief user's
#    guide (postscipt file guide.ps).

#######################################################################

#    test functions for examples in the Brief User's Guide   <<<<<<<<<<

#######################################################################


guide.example.tests.part1 <- function( verbose=T, synopsis=T, fuzz.small=1e-14, fuzz.large=1e-8, graphics=T, pause=F)
{# test examples in Brief User's guide
 # NOTE: it was necessary to reduce fuzz from 1e-14 because of differences
 # in the results between Splus 3.2 and Splus 3.3 (C libraries were changed).
 # Differences affected lsfit (used in est.VARX.ls) among other things.


  # If no device is active then write to postscript file 
  if (graphics)
   {if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
    else
      {old.par <- par()
       #on.exit(par(old.par))
      }
    if(pause) dev.ask(ask=T)
   }


  max.error <- NA
  if (synopsis & !verbose) cat("All Brief User Guide example part 1 tests ...")


  if (verbose) cat("Guide part 1 test 0 ... ")
    {if (is.R())
       {data(eg1.DSE.data.diff, package="dse1")
        data(eg1.DSE.data, package="dse1")
        data(egJofF.1dec93.data, package="dse1")
       } 
     if (is.S())
       {if(!exists("eg1.DSE.data.diff")) warning("eg1.DSE.data.diff does not exist")
        if(!exists("eg1.DSE.data"))      warning("eg1.DSE.data does not exist")
       } 
    }
  if (verbose) { cat("ok\n") }
     
  if (verbose) cat("Guide part 1 test 1 ... ")
  # previously search() was used to determine "from", but DSE.HOME is better
  #  if(1==pmatch("MS Windows",version$os, nomatch=0))
  #     from <- (search()[grep("b*/DSE/_Data", search())])[1]
  #  else
  #     from <- (search()[grep("b*/DSE/.Data", search())])[1]
  # if DSE is not in the search path then data file is assumed to be in 
  #  the present working directory
  #  if(5 < nchar(from)) from<-substring(from, first=1, last = nchar(from) -5)
  #  from <- paste(from,"eg1.dat", sep="")

  from <- paste(DSE.HOME, "/data/eg1.dat", sep="")
  data <- t(matrix(dsescan(from), 5,364))[,2:5]
  #  data <- list(
  #     input=tframed(data[,1  ,drop=F],  list(start=c(1961,3), frequency=12)),
  #     output=tframed(data[,2:4,drop=F], list(start=c(1961,3), frequency=12)))
  data <- TSdata(input=tframed(data[,1  ,drop=F],
                                         list(start=c(1961,3), frequency=12)),
              output=tframed(data[,2:4,drop=F], 
                                         list(start=c(1961,3), frequency=12)))
  input.series.names(data)   <-  "u1"
  output.series.names(data) <-  c("y1","y2","y3")
  error <- abs(126943980.50000011921 - sum(output.data(data)))
  ok <- 100*fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("Guide part 1 test 2 ... ")
  model1 <- est.VARX.ls(data, warn=F)
  model2 <- est.SS.Mittnik(data, n=14)
#  summary(model1)
#  summary(model2)
#  print(model1)
#  print(model2)
#  stability(model1)
#  stability(model2)
   if (graphics) tfplot(model1)

  error <- max(Mod(c(15.430979953081722655   - sum(TSmodel(model1)$A),
                         -1.1078692933906153506  - sum(TSmodel(model2)$F),
                         2.4561249653768193468   - sum(roots(model2)) )))
#                        2.4561249653768193468+0i- sum(roots(model2)) )))
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  ok <-  ok & is.TSestModel(model1) & is.TSestModel(model2)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("Guide part 1 test 3 ... ")
  ar<-array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2))
  ma<-array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2))
  arma<-ARMA(A=ar,B=ma,C=NULL)
#  print(arma)
  ok <- is.TSmodel(arma) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 4 ... ")
  data.arma.sim<-simulate(arma)
  arma<-l(arma,data.arma.sim)
#  summary(arma)
  if (graphics) 
     {tfplot(data.arma.sim)
      tfplot(arma)
     }
  ok <- is.TSdata(data) & is.TSestModel(arma) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 5 ... ")
  f<-array(c(.5,.3,.2,.4),c(2,2))
  h<-array(c(1,0,0,1),c(2,2))
  k<-array(c(.5,.3,.2,.4),c(2,2))
  ss<-SS(F=f,G=NULL,H=h,K=k)
#  ss
  ok <- is.SS(ss)
  all.ok <- all.ok & ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 6 ... ")
  data.ss.sim<-simulate(ss)
  ss<-l(ss,data.ss.sim)
#  summary(ss)
  if (graphics) tfplot(ss)
  ok <- is.TSestModel(ss)
  all.ok <- all.ok & ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("Guide part 1 test 7 ... ")
  ss.from.arma <- l(to.SS(arma), data.arma.sim)
  arma.from.ss <- l(to.ARMA(ss), data.ss.sim)
#  summary(ss.from.arma)
#  summary(arma)
#  summary(arma.from.ss)
#  summary(ss)
#  stability(arma)
#  stability(ss.from.arma)
#  caution: tests on $estimates will depend on seed when data is generated.
  error <- max(Mod(c(-0.15000000000000018874 - sum(TSmodel(ss.from.arma)$F),
                         0.47999999999999998224  - sum(TSmodel(arma.from.ss)$A),
                         -1                      - sum(roots(ss.from.arma)) )))
#                        -1+0i                   - sum(roots(ss.from.arma)) )))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (synopsis) 
    {if (verbose) cat("All Brief User Guide example part 1 tests completed")
     if (all.ok) cat(" OK\n")
     else    cat(", some FAILED! max.error = ", max.error,"\n")
    }
  invisible(all.ok)
}




guide.example.tests.part2 <- function( verbose=T, synopsis=T,
     fuzz.small=1e-14, fuzz.large=1e-8, graphics=T, pause=F, print.values=F)
{# test examples in Brief User's guide
  if (synopsis & !verbose) cat("Brief User Guide example part 2 tests ...")
  if (verbose) cat("Guide part 2 test 1 ... ")

  # If no device is active then write to postscript file 
  if (graphics)
   {if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
                {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
    else
      {old.par <- par()
       #on.exit(par(old.par))
      }
    if(pause) dev.ask(ask=T)
   }

  max.error <- NA
  all.ok <- T

  if (is.S())
    {test.rng1 <- list(kind="default", normal.kind="default", 
                       seed=c(13,44,1,25,56,0,6,33,22,13,13,0) )
     test.rng2 <- list(kind="default", normal.kind="default", 
                        seed=c(13,43,7,57,62,3,30,29,24,54,47,2) )
     test.rng4 <- list(kind="default", normal.kind="default", 
                        seed=c(29,55,47,18,33,1,15,15,34,46,13,2) )
     test.rng3 <- list(kind="default", normal.kind="default", 
                        seed=c( 53,41,26,39,10,1,19,25,56,32,28,3) )
    }
  else if (is.R()) 
    {test.rng1 <- test.rng2 <- test.rng3 <- test.rng4 <- 
           list(kind="Wichmann-Hill", normal.kind="Kinderman-Ramage",
	        seed=c(979, 1479, 1542))} # defaults changed in R 1.0.0
	   # might consider also 
	   #list(kind="Wichmann-Hill", normal.kind="user-supplied", seed=c(979, 1479, 1542))}
           # R's Box-Muller was declared not reproducible. 
	   
  from <- paste(DSE.HOME, "/data/eg1.dat", sep="")
  eg1.DSE.data <-example.get.eg.raw.data(from) #retrieves data from file
  eg2.DSE.data.names <- TSPADIdata(
	output="I37005", output.names="manuf.prod." , server="ets")
# Fame call disabled for testing:  eg2.DSE.data <- freeze(eg2.DSE.data.names)
  eg3.DSE.data.names <- TSPADIdata(
 	input="lfsa455", input.transforms= "percent.change",
 	input.names= "manuf.emp.",
 	output= "i37005", output.names=c( "manuf.prod."), 
        output.transforms= c("percent.change"),
 	  server="ets", pad.start=F, pad.end =T  )
# Fame call disabled for testing: eg3.DSE.data <- freeze(eg3.DSE.data.names)

 egJofF.1dec93.data.names <- TSPADIdata(
	input = "B14017", input.transforms= "diff", input.names="R90",
#	output = c("P484549", "I37026", "b1627", "b14013", discont.
	output = c("B820600", "I37026", "b1627", "b14013",
		   "b4237", "D767608", "b3400", "M.BCPI", "M.JQIND", "M.CUSA0"),
	output.transforms=c("percent.change", 
			"percent.change","percent.change",
			"diff", "diff", "percent.change",
			"percent.change", "percent.change",
			"percent.change", "percent.change"),
	output.names=c("CPI", "GDP", "M1", "RL", "TSE300", 
			"employment", "PFX", "com. price ind.", 
			"US ind. prod.", "US CPI"),
	server="ets")
# Fame call disabled for testing:
#   egJofF.1dec93.data <- freeze(egJofF.1dec93.data.names)
if      (is.R()) data("egJofF.1dec93.data", package="dse1")
else if (is.S()) source(paste(DSE.HOME, "/data/egJofF.1dec93.data.R", sep=""))

if(!exists("egJofF.1dec93.data"))warning("egJofF.1dec93.data does not exist")

  error <- abs(3352.4721630925987483 - 
            sum(c(output.data(egJofF.1dec93.data),
                   input.data(egJofF.1dec93.data))))
  ok <-  fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok &ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

# Section 4 - Model Estimation

  if (verbose) cat("Guide part 2 test 2 ... ")
  model.eg1.ls <- est.VARX.ls(trim.na(eg1.DSE.data), warn=F)
#  opts <-options(warn=-1) 
    subsample.data <- tfwindow(eg1.DSE.data,start=c(1972,1),end=c(1992,12),warn=F)
#  options(opts)
  # summary(model.eg1.ls)
  # print(model.eg1.ls)
  if (graphics)
    {tfplot(model.eg1.ls)
     tfplot(model.eg1.ls, start.=c(1990,1))
    }
  z <- check.residuals(model.eg1.ls, plot.=F, pac=T)
  if (is.S())      check.value <-
              c(4.67445135116577148, 3274.42578125,      -2371.9997808950302)
  else if (is.R()) check.value <- 
              c(4.674448837156188,   3274.422653488969,  -2371.999780895034)
# using my old acf instead of bats version gives
#             c(4.6744488371561879,      0.0,       -2371.999780895033837)

  if (print.values) print.test.value(c(sum(z$acf),sum(z$pacf),sum(z$cusum)) )  
  error <- max(abs(check.value   -   c(sum(z$acf),sum(z$pacf),sum(z$cusum)) ))
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 3 ... ")
  # NB- non-stationary data. ar is not really valid
  model.eg1.ar <- est.VARX.ar(trim.na(eg1.DSE.data), warn=F) 
  model.eg1.ss <- est.SS.from.VARX(trim.na(eg1.DSE.data), warn=F) 
# model.eg1.mle <- est.max.like(trim.na(eg1.DSE.data),model.eg1.ar) # this may be slow
  if (is.S())      check.value <- c(6738.642280883833, 6921.352391513382)
  else if (is.R()) check.value <- c(6738.562190977154, 6921.352391513382)#ts ar
  #elseif(is.R()) check.value <- c(6735.139112062216, 6921.352391513382)bats ar
# using my old ar:=ls gives      c(6921.352391513380, 6921.352391513380)#ar:=ls
  if (print.values) print.test.value(
      c(model.eg1.ar$estimates$like[1],model.eg1.ss$estimates$like[1]) )  
  error <- max(abs(check.value -
             c(model.eg1.ar$estimates$like[1],model.eg1.ss$estimates$like[1])))
  ok <- 10*fuzz.large > error    
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 4a... ")
   eg4.DSE.data<- egJofF.1dec93.data
   output.data(eg4.DSE.data) <- output.data(eg4.DSE.data, series=c(1,2,6,7))
 # following is optional 
 # tframe(output.data(eg4.DSE.data))<- tframe(output.data(egJofF.1dec93.data))

  model.eg4.bb <- est.black.box(trim.na(eg4.DSE.data), max.lag=3, verbose=F) 
  error <- abs(614.70500313590287078 - model.eg4.bb$estimates$like[1] )
  ok <-  fuzz.large > error 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 4b... ")
  z <- information.tests(model.eg1.ar, model.eg1.ss, Print=F, warn=F)
#  if (is.S())      check.value <- 231152.464267979725
#  else if (is.R()) check.value <- 231151.0300943982  # ts ar
# else if (is.R()) check.value <- 230978.2532793634  bats ar
# using my old ar=ls gives        233856.9237566061   #using ls for ar
# a small change in the accounting for degenerate subspaces in dse.2000.4 gives
  if (is.S())      check.value <- 225328.464267979754
  else if (is.R()) check.value <- 225327.03009431256
  if (print.values) print.test.value(sum(z[!is.na(z)]) )  
  error <- abs(check.value - sum(z[!is.na(z)]) )
  ok <-  1000*fuzz.large > error #fuzz.large works in Solaris but not Linux
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }


# Section 5 - Forecasting

  if (verbose) cat("Guide part 2 test 5 ... ")
  eg4.DSE.model <- est.VARX.ls(eg4.DSE.data)
#  Fame call disabled for testing: new.data <- freeze(eg4.DSE.data.names) 
  new.data <- TSdata(
              input= ts(rbind(input.data(eg4.DSE.data), matrix(.1,10,1)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)),    
              output=ts(rbind(output.data(eg4.DSE.data),matrix(.3,5,4)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)))
  series.names(new.data) <- series.names(eg4.DSE.data)
  z  <- l(TSmodel(eg4.DSE.model), trim.na(new.data)) 
#  z <- l(TSmodel(eg4.DSE.model), trim.na(freeze(eg4.DSE.data.names)))
  error <- max(abs(556.55870662521476788 -z$estimates$like[1]))
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 6 ... ")
  zz <- forecast(TSmodel(eg4.DSE.model), new.data)
  z <-  forecast(TSmodel(eg4.DSE.model), trim.na(new.data), 
		conditioning.inputs=input.data(new.data))
  if (graphics) tfplot(zz, start.=c(1990,6))
  error <- abs(4.7990339556773520258 - sum(forecasts(z)[[1]]))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  ok <- test.equal(zz,z) & ok
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 7 ... ")
  z <- forecast(eg4.DSE.model, conditioning.inputs.forecasts=matrix(.5,6,1)) 
  # Fame call disabled for testing: 
  # z <- forecast(TSmodel(eg4.DSE.model), freeze(egJofF.1dec93.data.names), 
  #		conditioning.inputs.forecasts=matrix(.5,6,1))
  # summary(z)
  # print(z)
  if (graphics)
    {tfplot(z)
     tfplot(z, start.=c(1990,1))
    }
  #  forecasts(z)[[1]]
  #  tfwindow(forecasts(z)[[1]], start=c(1994,5)) 
  ok <- all(start(eg4.DSE.model$data)   == c(1974,2)) &
        all(start(egJofF.1dec93.data) == c(1974,2)) 
  error <- max(abs(c(5.9414711908521793404 - sum(forecasts(z)[[1]][1:6,]),
                  3.7410224783909828972 - 
                     sum(tfwindow(forecasts(z)[[1]], start=c(1993,12), warn=F)))))
  ok <-  ok & (fuzz.small > error) 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 8 ... ")
  z <- l(TSmodel(eg4.DSE.model), new.data)
  if (graphics) tfplot(z)
  z <- feather.forecasts(TSmodel(eg4.DSE.model), new.data)
  if (graphics) tfplot(z)
  zz <-feather.forecasts(TSmodel(eg4.DSE.model), new.data,
                          from.periods =c(20,50,60,70,80), horizon=150)
  if (graphics) tfplot(zz)
  error <- max(abs(c(54.838475604100473504 -
                           sum( forecasts(z)[[1]][10:46,]),
                       53.824873541066445171 - 
                           sum(forecasts(zz)[[5]][80:116,]))))
  ok <-fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 9 ... ")
  z <- horizon.forecasts(TSmodel(eg4.DSE.model), new.data, horizons=c(1,3,6))
  if (graphics) tfplot(z)
#  error <- abs(653.329319170802592 - sum(z$horizon.forecasts) )
  error <- abs(653.329319170802592 - sum(forecasts(z)) )
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 10... ")
  fc1 <- forecast.cov(TSmodel(eg4.DSE.model), data=eg4.DSE.data)
  if (graphics) 
    {tfplot(fc1)
     tfplot(forecast.cov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, horizons= 1:4)) 
    }
  fc2 <- forecast.cov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, zero=T, trend=T)
  if (graphics) tfplot(fc2)
  error <- max(abs(c(14.933660144821400806 - sum(fc1$forecast.cov[[1]]),
                        14.933660144821400806 - sum(fc2$forecast.cov[[1]]),
                        31.654672476928297442 - sum(fc2$forecast.cov.zero),
                        18.324461923341953451 - sum(fc2$forecast.cov.trend) )))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  if (is.na(ok)) ok <- F
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 11... ")
  mod1 <- ARMA(A=array(c(1,-.25,-.05), c(3,1,1)), B=array(1,c(1,1,1)))
  mod2 <- ARMA(A=array(c(1,-.8, -.2 ), c(3,1,1)), B=array(1,c(1,1,1)))
  mod3 <- ARMA(
 	A=array(c( 
 	1.00,-0.06,0.15,-0.03,0.00,0.02,0.03,-0.02,0.00,-0.02,-0.03,-0.02,
	0.00,-0.07,-0.05,0.12,1.00,0.20,-0.03,-0.11,0.00,-0.07,-0.03,0.08,
 	0.00,-0.40,-0.05,-0.66,0.00,0.00,0.17,-0.18,1.00,-0.11,-0.24,-0.09 )
		,c(4,3,3)), 
 	B=array(diag(1,3),c(1,3,3)))
  e.ls.mod1 <- eval.estimation( mod1, replications=100, 
 	rng=test.rng1,
 	simulation.args=list(sampleT=100, sd=1), 
 	estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
 	criterion="TSmodel", quiet=T)

#    e.ar.mod1 <- eval.estimation( mod1, replications=100, 
#   	rng=test.rng1,
#   	simulation.args=list(sampleT=100, sd=1), 
#   	estimation="est.VARX.ar", estimation.args=list(max.lag=2, aic=F), 
#   	criterion="TSmodel", quiet=T)
#   tfplot(parms(e.ar.mod1))


  if (is.S())      check.value <- -0.29855874505752699744
  else if (is.R()) check.value <- -0.3699580622977686
  error <- abs( check.value - sum(parms(e.ls.mod1$result[[100]])))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 12... ")
  e.ls.mod2 <- eval.estimation( mod2, replications=100, 
                     rng=test.rng2,
                     simulation.args=list(sampleT=100, sd=1), 
                     estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
                     criterion="TSmodel", quiet=T)
  if (graphics)
    {old.par <- par(mfcol=c(2,1)) #set the number of plots on the plotics device
     on.exit(par(old.par))
     tfplot(parms(e.ls.mod1))
     tfplot(parms(e.ls.mod2)) 
     old.par <- c(old.par, par(mfcol=c(2,1)) )
     tfplot(parms(e.ls.mod1), cum=F, bounds=F) 
     tfplot(parms(e.ls.mod2), cum=F, bounds=F) 
     distribution(parms(e.ls.mod1), bandwidth=.2)
     distribution(parms(e.ls.mod2), bandwidth=.2)
    } 

  if (is.S())      check.value <- -1.0021490287427212706
  else if (is.R()) check.value <- -1.0028944627996934
  error <- abs(check.value - sum(parms(e.ls.mod2$result[[100]])))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 13... ")
  e.ls.mod1.roots <- roots(e.ls.mod1)
  if (graphics)
    {plot(e.ls.mod1.roots) 
     plot(e.ls.mod1.roots, complex.plane=F)
     plot(roots(e.ls.mod2), complex.plane=F) 
     distribution(e.ls.mod1.roots, bandwidth=.2) 
     distribution(roots(e.ls.mod2), bandwidth=.1) 
    }

  if (is.S())      check.value <- 0.36159459310761993267
  else if (is.R()) check.value <- 0.2119677206564640
# error <- Mod(0.36159459310761993267+0i - sum(e.ls.mod1.roots$result[[100]]))
  error <- Mod(check.value   - sum(e.ls.mod1.roots$result[[100]]))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 14... ")
  pc <- forecast.cov.estimators.wrt.true(mod3,
 	rng=test.rng3,
 	estimation.methods=list(est.VARX.ls=list(max.lag=6)),
 	est.replications=2, pred.replications=10, quiet=T)
  # the fuzz.small has to be relaxed here to accomodate differences in rnorm
  #   between Splus3.1 and Splus3.2  (the numbers are from Splus3.2)

  if      (is.S())  check.value <-
      c(60.927013860328429473, 62.32729288591478678, 63.17808145947956433) 
  else if (is.R()) check.value <-
      c( 54.164759056117504, 53.519297277839669, 59.341526159411558)
  error <- max(abs(check.value -c(sum(pc$forecast.cov[[1]]), 
                      sum(pc$forecast.cov.zero), sum(pc$forecast.cov.trend) )))
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 15... ")
  pc.rd <- forecast.cov.reductions.wrt.true(mod3,
 	rng=test.rng4,
 	estimation.methods=list(est.VARX.ls=list(max.lag=3)),
 	est.replications=2, pred.replications=10, quiet=T)

  if (is.S())      check.value <-
         c(58.75543799264762157,60.451513998215133938, 64.089618782185240775) 
  else if (is.R()) check.value <- 
         c( 51.237201863944890, 53.519297277839669, 59.341526159411558)
  error <- max(abs(check.value - c(sum(pc.rd$forecast.cov[[1]]),
                  sum(pc.rd$forecast.cov.zero), sum(pc.rd$forecast.cov.trend))))
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 2 test 16... ")
  z <-out.of.sample.forecast.cov.estimators.wrt.data(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.VARX.ar=list(warn=F), 
 		est.VARX.ls=list(warn=F)), 
 	trend=T, zero=T)
  if (graphics) tfplot(z)
  opts <- options(warn=-1)
  zz <-out.of.sample.forecast.cov.estimators.wrt.data(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.black.box4=list(max.lag=3, verbose=F, warn=F),
		est.VARX.ls=list(max.lag=3, warn=F)), 
	trend=T, zero=T)

#    zf<-horizon.forecasts(zz$multi.model[[1]],zz$data, horizons=c(1,3,6))
    zf<-horizon.forecasts(TSmodel(zz, select=1),TSdata(zz), horizons=c(1,3,6))
  options(opts)
  zf<- zf$horizon.forecasts[3,30,]
  if (graphics) tfplot(z)
  if (is.S())      check.value <- 
     c(6120.97621905043979, 175568.040899036743, 24.568074094041549,
       1e-10*c(158049871127.845642, 3330592793.50789356, 
              1242727188.69001055, 1606263575.00784183))

  else if (is.R()) 
    {check.value <-
     c(6120.97621905044,  175568.0408990367,  24.56807409404155, 
       1e-10*c(158034797997.4372,  3330592793.507894,
              1242727188.690011,  1606263575.007842))
# using my old ls for ar instead of bats version gives
#    c(6120.9762190509673, 175568.04089903546, 24.568074093999403,
#       1e-10*c(3330592793.507894, 3330592793.507894,
#              1242727188.690011, 1606263575.0078418)) # using ls for ar
    }

  if (print.values) print.test.value(c(zf,  # 1e-5*
                     c(sum( z$forecast.cov[[1]]), sum( z$forecast.cov[[2]]),
                       sum(zz$forecast.cov[[1]]), sum(zz$forecast.cov[[2]]))) )  
  error <- max(abs(check.value -
         c(zf,  1e-10*c(sum( z$forecast.cov[[1]]), sum( z$forecast.cov[[2]]),
                       sum(zz$forecast.cov[[1]]), sum(zz$forecast.cov[[2]])))))
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (synopsis) 
    {if (verbose) cat("All Brief User Guide example tests part 2 completed")
     if (all.ok) cat(" OK\n")
     else    cat(", some FAILED! max.error = ", max.error,"\n")
    }
  invisible(all.ok)
}


guide.example.tests.part3 <- function( verbose=T, synopsis=T, fuzz.small=1e-14, fuzz.large=1e-8, graphics=T, pause=F)
{
 if (!exists("check.padi.server"))
    stop("This test requires a PADI server called ets and is primarily intended for checks at the BOC. The code can be examined as an example.")
 else  if (!check.padi.server("ets")) # ets bc
    stop("This test requires a server called ets and is primarily intended for checks at the BOC. The code can be examined as an example.")

  # If no device is active then write to postscript file 
     if (graphics)
   {if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
                {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
    else
      {old.par <- par()
       #on.exit(par(old.par))
      }
    if(pause) dev.ask(ask=T)
   }

  max.error <- NA
  if (synopsis & !verbose) cat("Brief User Guide example part 3 tests ...")
  if (verbose) cat("Guide part 3 test 1 ... ")

# help.start.DSE(browser="mosaic")
# help.start.DSE()
# openlook()  #for plot window


# manufacturing

#  two outputs, no input

   cbps.manuf.data2.ids <- TSPADIdata(
      output=c("i37013", "i37005"), output.transforms="percent.change",
      output.names=c("cbps.prod.", "manuf.prod."),
       server="ets", pad.start=F, pad.end =T )
   # get data
   cbps.manuf.data2 <- freeze(cbps.manuf.data2.ids)


   # Estimate models with employment as input (exogenous) variable.

   manuf.data.ids <- TSPADIdata(
      input ="lfsa455", input.transforms="percent.change",
      input.names="manuf.emp.",
      output="I37005", output.transforms="percent.change",
      output.names="manuf.prod.",
      server="ets", pad.start=F, pad.end =T  )
      #pad.end=F for estimation or use trim.na(data) 
   manuf.data <- freeze(manuf.data.ids)
   if (graphics) tfplot(manuf.data)       # shows NA data in the middle
   manuf.data <- tfwindow(manuf.data, start=c(1976,2))
   manuf.data.ids <- modify.TSPADIdata(manuf.data.ids,
                     start=c(1976,2)) # avoid NAs in middle as with window above
   manuf.data <- freeze(manuf.data.ids)

   #if (graphics) tfplot(tfwindow(manuf.data, start=c(1995,11)))   
   if (graphics) tfplot(manuf.data, start.=c(1995,11))  # Bug the . after start is nec.
   # cbps
   #  two inputs
   cbps.manuf.data.ids <- TSPADIdata(
      input =c("lfsa462","lfsa455"),
      input.transforms="percent.change", 
      input.names=c("cbps.emp.", "manuf.emp"),
      output="i37013",
      output.transforms="percent.change",
      output.names="cbps.prod.",
      start=c(1976,2),
      server="ets", db="",
      pad.start=F,  pad.end =T  )
#   cbps.manuf.data.ids <- modify(cbps.manuf.data.ids, start=c(1976,2) ) # avoid NAs in middle
   cbps.manuf.data <- freeze(cbps.manuf.data.ids)
   # cbps.manuf.data <- tfwindow(cbps.manuf.data, start=c(1976,2))

#  two outputs, one input

   cbps.manuf.data3.ids <- TSPADIdata(
      input ="lfsa462",
      input.transforms="percent.change",input.names="cbps.emp.",
      output=c("i37013", "i37005"), 
      output.transforms=c("percent.change", "percent.change"), 
      output.names=c("cbps.prod.","manuf.prod."),
      start=c(1976,2),
      server="ets", db ="", pad.start=F, pad.end =T)
#   cbps.manuf.data3.ids <- modify(cbps.manuf.data3.ids, start=c(1976,2)) 
   cbps.manuf.data3 <- freeze(cbps.manuf.data3.ids)
   # cbps.manuf.data3 <- tfwindow(cbps.manuf.data3, start=c(1976,2))

  ok <- is.TSdata(manuf.data) & is.TSdata(cbps.manuf.data) & 
        is.TSdata(cbps.manuf.data3 )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 3 test 2 ... ")

# estimate model

# dev.ask(T)

# manuf.model  <- bft(trim.na(manuf.data)) # might take some time
# or
# manuf.model  <- bft(trim.na(manuf.data), max.lag=5) 
# or
manuf.model  <- bft(trim.na(manuf.data), verbose=F, max.lag=5)
# may give "Warning: Cannot open audit file"  (this is not a problem)


 #   manuf.model  # display the model parameters

   if (graphics) tfplot(manuf.model)
   if (graphics) tfplot(manuf.model, start.=c(1990,1))
   if (graphics) tfplot(manuf.model, start.=c(1995,1))
   cbps.manuf.model <- bft(trim.na(cbps.manuf.data),verbose=F)
   if (graphics) tfplot(cbps.manuf.model)
   if (graphics) tfplot(cbps.manuf.model, start.=c(1995,1))

# to forecast with the model using all available employment data
   zd <- trim.na(manuf.data)
   if (all(end(zd) == end(trim.na(input.data(manuf.data)))))
     zd <- tfwindow(zd, end=end(zd) -c(0,1)) # otherwise no forecast is produced
   z <- forecast(TSmodel(manuf.model), zd,   
          conditioning.inputs=trim.na(input.data(manuf.data)))
   if (graphics) tfplot(z, start.=c(1995,1))

# to see the forecast

   zot <- forecasts(z)[[1]]
   zot <- tfwindow(forecasts(z)[[1]], start=c(1996,3), warn=F)


  ok <- is.forecast(z)
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  if (verbose) cat("Guide part 3 test 3 ... ")

   fc <- forecast.cov(manuf.model)
   if (graphics) tfplot(fc)                                

#relative to zero and trend ...

# to analyse forecast errors

   fc <- forecast.cov(manuf.model,  zero=T, trend=T)
   if (graphics) tfplot(fc)   

# analysing out-of-sample forecasts errors

#   outfc <-out.of.sample.forecast.cov.estimators.wrt.data(trim.na(manuf.data),
#    estimation.sample=.5,
#    estimation.methods = list(bft=list(verbose=F)), trend=T, zero=T)

   outfc <-out.of.sample.forecast.cov.estimators.wrt.data(trim.na(manuf.data),
    estimation.sample=.5,
    estimation.methods = list(bft=list(verbose=F), est.VARX.ls=NULL),
    trend=T, zero=T)

#   outfc <-out.of.sample.forecast.cov.estimators.wrt.data(
#    trim.na(cbps.manuf.data3),
#    estimation.sample=.5,
#    estimation.methods = list(bft=list(verbose=F), est.VARX.ls=NULL),
#    trend=T, zero=T)

   if (graphics) tfplot(outfc)



# Sometimes it is useful to send output to a text file rather than to the 
#  screen. To do this use

   sink(file= "zzz.some.name")  #all output goes to file
   ls()
   sink()    # output returns to normal
   syskern.rm("zzz.some.name")

   new.data <- freeze(manuf.data.ids)  # retrieve new (updated) data set

# To run the model with different data you can use

   z <- l(TSmodel(manuf.model), trim.na(new.data)) 
   z <- l(TSmodel(manuf.model), trim.na(freeze(manuf.data.ids)))
   if (graphics) tfplot(z)
   if (graphics) tfplot(z, start.=c(1995,8))
#   z <- forecast(TSmodel(manuf.model), trim.na(new.data),  
#     conditioning.inputs=input.data(new.data))
# These series seem to have been discontinued in 1996,12 and 1997,9 so the 
#   following trunction is done for testing
   z <- forecast(TSmodel(manuf.model), trim.na(tfwindow(new.data, end=c(1996,8))), 
     conditioning.inputs=input.data(tfwindow(new.data, end=c(1996,12))))

   if (graphics) tfplot(z, start.=c(1995,6))

#if you actually want the numbers type

   zot <- forecasts(z)[[1]]

   zot <- tfwindow(forecasts(z)[[1]], start=c(1996,2), warn=F) 
# to put the projected data into a Fame database (at BOC using PADI)

   zot <- output.series.names(z)

   putpadi(forecasts(z)[[1]],dbname="zzz.nameofdatabase.db", 
        series=output.series.names(z))
   syskern.rm("zzz.nameofdatabase.db")

  ok <- T
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }


  if (verbose) cat("Guide part 3 test 4 ... ")

# To give the initial set up before the first time this is run

     manuf.previous.data <- freeze(manuf.data.ids)
  # then to test (so it looks like the data has changed):
     manuf.previous.data$output[1,1] <- NA

warning("Skipping part 3 test 4, data has been discontinued.")
#   r <-simple.monitoring(manuf.model, 
#           manuf.data.ids, 
#           manuf.previous.data,
#           mail.list=user.name(),
#           message.title="    Manufacturing Monitoring TEST !!!! ",
#           message.subject=" TEST !!!! Manufacturing Monitoring",
#           show.start= c(0,-3), 
#           report.variables=series.names(manuf.data.ids),
#           data.sub.heading=
#               "   %chg       %chg",
#           message.footnote="              f - forecast value" ,
#           data.tag=" ",
#           forecast.tag="f" )
        # ,save.as=paste("Manufacturing.monitoring.",
        #               make.names(paste( date.parsed(), collapse=".")),sep="")
        # ,run.again=F
#   manuf.previous.data <- r[["data"]]
#   zot <- r$status

  ok <- T
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }


  if (synopsis) 
    {if (verbose) cat("All Brief User Guide example tests part 3 completed")
     if (all.ok) cat(" OK\n")
     else    cat(", some FAILED! max.error = ", max.error,"\n")
    }
  invisible(all.ok)
}


#######################################################################

#                    end

#######################################################################



#######################################################################

#    test functions for examples in the Brief User's Guide   <<<<<<<<<<

#######################################################################



# The tests requiring padi and dsepadi are in dsepadi/tests.

require("dse2")
 Sys.info()
 version.dse()

print.values <- F
fuzz.small <- 1e-14
fuzz.large <- 1e-8
max.error <- NA
all.ok <- T


  cat("Guide part 2 test 1 ... \n")

  if (is.S())
    {test.rng1 <- list(kind="default", normal.kind="default", 
                       seed=c(13,44,1,25,56,0,6,33,22,13,13,0) )
     test.rng2 <- list(kind="default", normal.kind="default", 
                        seed=c(13,43,7,57,62,3,30,29,24,54,47,2) )
     test.rng4 <- list(kind="default", normal.kind="default", 
                        seed=c(29,55,47,18,33,1,15,15,34,46,13,2) )
     test.rng3 <- list(kind="default", normal.kind="default", 
                        seed=c( 53,41,26,39,10,1,19,25,56,32,28,3) )
    } else if (is.R()) 
    {test.rng1 <- test.rng2 <- test.rng3 <- test.rng4 <- 
           list(kind="Wichmann-Hill", normal.kind="Kinderman-Ramage",
	        seed=c(979, 1479, 1542))} # defaults changed in R 1.0.0
	   # might consider also 
	   #list(kind="Wichmann-Hill", normal.kind="user-supplied", seed=c(979, 1479, 1542))}
           # R's Box-Muller was declared not reproducible. 
	   

eg1.DSE.data <- t(matrix(scan(paste(DSE.HOME, "/data/eg1.dat", sep="")),
                         5, 364))[, 2:5] 

eg1.DSE.data <- TSdata(input = eg1.DSE.data[,1,drop = F], 
                      output = eg1.DSE.data[, 2:4, drop = F])
		      
eg1.DSE.data <-tframed(eg1.DSE.data, list(start=c(1961,3), frequency=12))

seriesNamesInput(eg1.DSE.data) <- "R90"
seriesNamesOutput(eg1.DSE.data) <- c("M1","GDPl2", "CPI")


if (is.R()) data("egJofF.1dec93.data", package="dse1")
if (is.S()) 
   {source(paste(DSE.HOME, "/data/egJofF.1dec93.data.R", sep=""))
    class(egJofF.1dec93.data$output) <- class(egJofF.1dec93.data$input) <- NULL
    }

if(!exists("egJofF.1dec93.data"))warning("egJofF.1dec93.data does not exist")

  error <- abs(3352.4721630925987483 - 
            sum(c(output.data(egJofF.1dec93.data),
                   input.data(egJofF.1dec93.data))))
  ok <-  fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok &ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

# Section 4 - Model Estimation

  cat("Guide part 2 test 2 ... \n")
  model.eg1.ls <- est.VARX.ls(trim.na(eg1.DSE.data), warn=F)
#  opts <-options(warn=-1) 
    subsample.data <- tfwindow(eg1.DSE.data,start=c(1972,1),end=c(1992,12),warn=F)
#  options(opts)
  # summary(model.eg1.ls)
  # print(model.eg1.ls)
  
  tfplot(model.eg1.ls)
  tfplot(model.eg1.ls, start.=c(1990,1))

 
  z <- check.residuals(model.eg1.ls, plot.=F, pac=T)
  check.value <- (
    if (is.S())      
              c(4.67445135116577148, 3274.42578125,      -2371.9997808950302)
  else if (is.R()) 
              c(4.674448837156188,   3274.422653488969,  -2371.999780895034) )
# using my old acf instead of bats version gives
#             c(4.6744488371561879,      0.0,       -2371.999780895033837)

  tst <- c(sum(z$acf),sum(z$pacf),sum(z$cusum))
  print.test.value(c(tst), digits=18)
  error <- max(abs(check.value   -   tst ))

  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 3 ... \n")
  # NB- non-stationary data. ar is not really valid
  model.eg1.ar <- est.VARX.ar(trim.na(eg1.DSE.data), warn=F) 
  model.eg1.ss <- est.SS.from.VARX(trim.na(eg1.DSE.data), warn=F) 
# model.eg1.mle <- est.max.like(trim.na(eg1.DSE.data),model.eg1.ar) # this may be slow
  if (is.S()) check.value <- c(6738.642280883833, 6921.352391513382) 
  if (is.R()) check.value <- c(6738.562190977154, 6921.352391513382)#ts ar
  #elseif(is.R()) check.value <- c(6735.139112062216, 6921.352391513382)bats ar
# using my old ar:=ls gives      c(6921.352391513380, 6921.352391513380)#ar:=ls

  tst <- c(model.eg1.ar$estimates$like[1], model.eg1.ss$estimates$like[1])  
  print.test.value(c(tst), digits=18)
  error <- max(abs(check.value - tst))

  ok <- 10*fuzz.large > error 
  if ((Sys.info()[["sysname"]] == "Linux") && ! ok) {
    warning("Using relaxed tolerance for Linux.")   
    ok <- 100*fuzz.large > error
    }   
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 4a... \n")
   eg4.DSE.data<- egJofF.1dec93.data
   output.data(eg4.DSE.data) <- output.data(eg4.DSE.data, series=c(1,2,6,7))
 # following is optional 
 # tframe(output.data(eg4.DSE.data))<- tframe(output.data(egJofF.1dec93.data))

  model.eg4.bb <- est.black.box(trim.na(eg4.DSE.data), max.lag=3, verbose=F) 

  tst <- model.eg4.bb$estimates$like[1]
  print.test.value(c(tst), digits=18)
  error <- abs(614.70500313590287078 - tst )

  ok <-  fuzz.large > error 
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }


  cat("Guide part 2 test 4b... \n")
  z <- information.tests(model.eg1.ar, model.eg1.ss, Print=F, warn=F)
#  if (is.S())      check.value <- 231152.464267979725
#  else if (is.R()) check.value <- 231151.0300943982  # ts ar
# else if (is.R()) check.value <- 230978.2532793634  bats ar
# using my old ar=ls gives        233856.9237566061   #using ls for ar
# a small change in the accounting for degenerate subspaces in dse.2000.4 gives
  if (is.S())      check.value <- 225328.464267979754 else
       if (is.R()) check.value <- 225327.03009431256
  if (print.values) print.test.value(sum(z[!is.na(z)]) )  
  error <- abs(check.value - sum(z[!is.na(z)]) )
  ok <-  1000*fuzz.large > error #fuzz.large works in Solaris but not Linux
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }


# Section 5 - Forecasting

  cat("Guide part 2 test 5 ... \n")
  eg4.DSE.model <- est.VARX.ls(eg4.DSE.data)
#  Fame call disabled for testing: new.data <- freeze(eg4.DSE.data.names) 
  new.data <- TSdata(
              input= ts(rbind(input.data(eg4.DSE.data), matrix(.1,10,1)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)),    
              output=ts(rbind(output.data(eg4.DSE.data),matrix(.3,5,4)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)))
  seriesNames(new.data) <- seriesNames(eg4.DSE.data)
  z  <- l(TSmodel(eg4.DSE.model), trim.na(new.data)) 
#  z <- l(TSmodel(eg4.DSE.model), trim.na(freeze(eg4.DSE.data.names)))
  error <- max(abs(556.55870662521476788 -z$estimates$like[1]))
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 6 ... \n")
  zz <- forecast(TSmodel(eg4.DSE.model), new.data)
  z <-  forecast(TSmodel(eg4.DSE.model), trim.na(new.data), 
		conditioning.inputs=input.data(new.data))
  tfplot(zz, start.=c(1990,6))
  error <- abs(4.7990339556773520258 - sum(forecasts(z)[[1]]))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  ok <- test.equal(zz,z) & ok
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 7 ... \n")
  z <- forecast(eg4.DSE.model, conditioning.inputs.forecasts=matrix(.5,6,1)) 
  # Fame call disabled for testing: 
  # z <- forecast(TSmodel(eg4.DSE.model), freeze(egJofF.1dec93.data.names), 
  #		conditioning.inputs.forecasts=matrix(.5,6,1))
  # summary(z)
  # print(z)
  tfplot(z)
  tfplot(z, start.=c(1990,1))
  
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
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 8 ... \n")
  z <- l(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  z <- featherForecasts(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  zz <-featherForecasts(TSmodel(eg4.DSE.model), new.data,
                          from.periods =c(20,50,60,70,80), horizon=150)
  tfplot(zz)
  error <- max(abs(c(54.838475604100473504 -
                           sum( forecasts(z)[[1]][10:46,]),
                       53.824873541066445171 - 
                           sum(forecasts(zz)[[5]][80:116,]))))
  ok <-fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 9 ... \n")
  z <- horizonForecasts(TSmodel(eg4.DSE.model), new.data, horizons=c(1,3,6))
  tfplot(z)
#  error <- abs(653.329319170802592 - sum(z$horizonForecasts) )
  error <- abs(653.329319170802592 - sum(forecasts(z)) )
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 10... \n")
  fc1 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data)
  
  tfplot(fc1)
  tfplot(forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, horizons= 1:4)) 
 
  fc2 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, zero=T, trend=T)
  tfplot(fc2)
  error <- max(abs(c(14.933660144821400806 - sum(fc1$forecastCov[[1]]),
                        14.933660144821400806 - sum(fc2$forecastCov[[1]]),
                        31.654672476928297442 - sum(fc2$forecastCov.zero),
                        18.324461923341953451 - sum(fc2$forecastCov.trend) )))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  if (is.na(ok)) ok <- F
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 11... \n")
  mod1 <- ARMA(A=array(c(1,-.25,-.05), c(3,1,1)), B=array(1,c(1,1,1)))
  mod2 <- ARMA(A=array(c(1,-.8, -.2 ), c(3,1,1)), B=array(1,c(1,1,1)))
  mod3 <- ARMA(
 	A=array(c( 
 	1.00,-0.06,0.15,-0.03,0.00,0.02,0.03,-0.02,0.00,-0.02,-0.03,-0.02,
	0.00,-0.07,-0.05,0.12,1.00,0.20,-0.03,-0.11,0.00,-0.07,-0.03,0.08,
 	0.00,-0.40,-0.05,-0.66,0.00,0.00,0.17,-0.18,1.00,-0.11,-0.24,-0.09 )
		,c(4,3,3)), 
 	B=array(diag(1,3),c(1,3,3)))
  e.ls.mod1 <- EstEval( mod1, replications=100, 
 	rng=test.rng1,
 	simulation.args=list(sampleT=100, sd=1), 
 	estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
 	criterion="TSmodel", quiet=T)

#    e.ar.mod1 <- EstEval( mod1, replications=100, 
#   	rng=test.rng1,
#   	simulation.args=list(sampleT=100, sd=1), 
#   	estimation="est.VARX.ar", estimation.args=list(max.lag=2, aic=F), 
#   	criterion="TSmodel", quiet=T)
#   tfplot(coef(e.ar.mod1))


  if (is.S()) check.value <- -0.29855874505752699744
  if (is.R()) check.value <- -0.3699580622977686
  error <- abs( check.value - sum(coef(e.ls.mod1$result[[100]])))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 12... \n")
  e.ls.mod2 <- EstEval( mod2, replications=100, 
                     rng=test.rng2,
                     simulation.args=list(sampleT=100, sd=1), 
                     estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
                     criterion="TSmodel", quiet=T)
     old.par <- par(mfcol=c(2,1)) #set the number of plots on the plotics device
     on.exit(par(old.par))
     tfplot(coef(e.ls.mod1))
     tfplot(coef(e.ls.mod2)) 
     old.par <- c(old.par, par(mfcol=c(2,1)) )
     tfplot(coef(e.ls.mod1), cum=F, bounds=F) 
     tfplot(coef(e.ls.mod2), cum=F, bounds=F) 
     distribution(coef(e.ls.mod1), bandwidth=.2)
     distribution(coef(e.ls.mod2), bandwidth=.2)
    

  if (is.S()) check.value <- -1.0021490287427212706
  if (is.R()) check.value <- -1.0028944627996934
  error <- abs(check.value - sum(coef(e.ls.mod2$result[[100]])))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 13... \n")
  e.ls.mod1.roots <- roots(e.ls.mod1)
  
     plot(e.ls.mod1.roots) 
     plot(e.ls.mod1.roots, complex.plane=F)
     plot(roots(e.ls.mod2), complex.plane=F) 
     distribution(e.ls.mod1.roots, bandwidth=.2) 
     distribution(roots(e.ls.mod2), bandwidth=.1) 
  

  if (is.S()) check.value <- 0.36159459310761993267
  if (is.R()) check.value <- 0.2119677206564640
# error <- Mod(0.36159459310761993267+0i - sum(e.ls.mod1.roots$result[[100]]))
  error <- Mod(check.value   - sum(e.ls.mod1.roots$result[[100]]))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 14... \n")
  pc <- forecastCov.estimatorsWRTtrue(mod3,
 	rng=test.rng3,
 	estimation.methods=list(est.VARX.ls=list(max.lag=6)),
 	est.replications=2, pred.replications=10, quiet=T)
  # the fuzz.small has to be relaxed here to accomodate differences in rnorm
  #   between Splus3.1 and Splus3.2  (the numbers are from Splus3.2)

  if (is.S())  check.value <-
      c(60.927013860328429473, 62.32729288591478678, 63.17808145947956433) 
  if (is.R()) check.value <-
      c( 54.164759056117504, 53.519297277839669, 59.341526159411558)
  error <- max(abs(check.value -c(sum(pc$forecastCov[[1]]), 
                      sum(pc$forecastCov.zero), sum(pc$forecastCov.trend) )))
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 15... \n")
  pc.rd <- forecastCov.reductionsWRTtrue(mod3,
 	rng=test.rng4,
 	estimation.methods=list(est.VARX.ls=list(max.lag=3)),
 	est.replications=2, pred.replications=10, quiet=T)

  if (is.S()) check.value <-
         c(58.75543799264762157,60.451513998215133938, 64.089618782185240775) 
  if (is.R()) check.value <- 
         c( 51.237201863944890, 53.519297277839669, 59.341526159411558)
  error <- max(abs(check.value - c(sum(pc.rd$forecastCov[[1]]),
                  sum(pc.rd$forecastCov.zero), sum(pc.rd$forecastCov.trend))))
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("Guide part 2 test 16... \n")
  z <-out.of.sample.forecastCov.estimatorsWRTdata(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.VARX.ar=list(warn=F), 
 		est.VARX.ls=list(warn=F)), 
 	trend=T, zero=T)
  tfplot(z)
  opts <- options(warn=-1)
  zz <-out.of.sample.forecastCov.estimatorsWRTdata(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.black.box4=list(max.lag=3, verbose=F, warn=F),
		est.VARX.ls=list(max.lag=3, warn=F)), 
	trend=T, zero=T)

#    zf<-horizonForecasts(zz$multi.model[[1]],zz$data, horizons=c(1,3,6))
    zf<-horizonForecasts(TSmodel(zz, select=1),TSdata(zz), horizons=c(1,3,6))
  options(opts)
  zf<- zf$horizonForecasts[3,30,]
  tfplot(z)
  if (is.S())      check.value <- 
     c(6120.97621905043979, 175568.040899036743, 24.568074094041549,
       1e-10*c(158049871127.845642, 3330592793.50789356, 
              1242727188.69001055, 1606263575.00784183))

  if (is.R()) 
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
                     c(sum( z$forecastCov[[1]]), sum( z$forecastCov[[2]]),
                       sum(zz$forecastCov[[1]]), sum(zz$forecastCov[[2]]))) )  
  error <- max(abs(check.value -
         c(zf,  1e-10*c(sum( z$forecastCov[[1]]), sum( z$forecastCov[[2]]),
                       sum(zz$forecastCov[[1]]), sum(zz$forecastCov[[2]])))))
  ok <-  fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  {if (ok) cat("ok\n") else cat("failed! error= ", error,"\n") }

  cat("All Brief User Guide example tests part 2 completed")
     if (all.ok) cat(" OK\n")  else 
        cat(", some FAILED! max.error = ", max.error,"\n")

  if (!all.ok) stop("Some tests FAILED")





# retrieve data from file eg1.dat
# Define some example generic model evaluation functions and
# define data verification functions and
# a suite of test functions (for DSE functions)
# These can be executed by
#    example.verify.data()
#    example.verify.raw.data()
#    example.tests()

example.get.eg.raw.data <- function(file)
{ example.raw.data <- t(matrix(dsescan(file), 5, 364))[, 2:5]
  example.raw.data <-list(
    input= ts(example.raw.data[,1,drop = F],
       start=c(1961,3), frequency=12),
    output=ts(example.raw.data[, 2:4, drop = F],
       start=c(1961,3), frequency=12) )
  dimnames(example.raw.data$input) <- list(NULL, "R90")
  dimnames(example.raw.data$output) <- list(NULL, c("M1","GDPl2", "CPI"))
  TSdata(example.raw.data)
}

example.convert.eg.data <- function(example.raw.data)
{ example.data <-list(
    input =ts(input.data(example.raw.data)[2:364, , drop = F], 
       start=c(1961,3), frequency=12),
    output=ts(diff(log(output.data(example.raw.data)[,, drop = F])),
       start=c(1961,3), frequency=12) )
  dimnames(example.data$input) <- list(NULL, "R90")
  dimnames(example.data$output) <- list(NULL, c("M1","GDPl2", "CPI"))
  TSdata(example.data)
}

example.verify.data <- function(example.data, verbose=T, fuzz.small=1e-14)
{# verify that data in example.data is correct data.
  z <-cbind(input.data(example.data),output.data(example.data))
  if (verbose) cat("example.data\n sample size...")
  ok <- 363 == dim(z)[1]
  ok <- ok & ( 4 == dim(z)[2] )
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }
  if (verbose) cat(" sample mean...")
  ok <- ok & all(fuzz.small > abs(apply(z, 2, mean) -
           c(8.606143250688707, 0.005502904131200907, 0.003297288176061855,
              0.004576853034062842)))
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }

  if (verbose) cat(" sample var ...")
  ok <- ok & all(fuzz.small > abs(apply(z, 2, var) -
             c(12.5442812169916, 0.0001384487711077435, 3.572716474599408e-05,
                 1.396066119144909e-05)))
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }
invisible(ok)
}

example.verify.raw.data <- function(example.raw.data, verbose=T, fuzz.small=1e-5)
{# verify that data in example.raw.data is correct data.
  z <-cbind(input.data(example.raw.data),output.data(example.raw.data))
  if (verbose) cat("example.raw.data\n sample size...")
  ok <- 364 == dim(z)[1]
  ok <- ok & ( 4 == dim(z)[2] )
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }
  if (verbose) cat(" sample mean...")
  ok <- ok & all(fuzz.small > abs(apply(z, 2, mean) -c(8.592884615384618, 19217.14560439561, 329471.3387362638, 58.71483516483516)))
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }

  if (verbose) cat(" sample var ...")
  ok <- ok & all(fuzz.small > abs(apply(z, 2, var) -c(12.57371204174613, 125360567.0779145, 11249376720.93049, 1067.923691157328)))
  if (verbose) 
    {if (ok) cat("ok\n")
     else cat("NOT CORRECT!\n")
    }
invisible(ok)
}

example.like.sub.sample <- function(model)
{#  example of generic evaluation using class and methods approach.
 # this example function takes a TSmodel and evaluates it using a 
 # fixed data set (the first 240 observations from example.data).
 l(model, example.data, sampleT=240, predictT=240)
}

example.like.total.sample <- function(model)
{#  example of generic evaluation using class and methods approach.
 # this example function takes a TSmodel and evaluates it using a 
 # fixed data set (example.data).
 l(model, example.data, sampleT=363, predictT=363)
}

example.tests <- function(example.data, verbose=T, summary=T,
        fuzz.small=1e-14, fuzz.large=1e-8)
{# A short set of tests of DSE functions using example.data.
 # Use as follows:
 #  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
 #  else if (is.S()) source(paste(DSE.HOME,"/data/eg1.DSE.data.diff.R", sep=""))
 #   example.tests(eg1.DSE.data.diff) 
 # The main short coming of these tests is that they do not test
 # functions which produce output, such as display, summary, graph
 # and check.residuals.
 # Note- The first part of these tests uses the total sample. 
 #    Working paper estimated with sub-sample!

  if (!example.verify.data(example.data,fuzz.small=fuzz.small, verbose=verbose))
     stop("example.data does not verify correctly. Example testing stopped.")
  max.error <- NA
  if (verbose) cat("Testing DSE functions with example.data.\n")
  if (verbose) cat(" est.VARX.ar...")
  VARmodel <- est.VARX.ar(example.data, re.add.means=F)
  error <- abs(-3879.7321062329338019 - VARmodel$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat(" est.VARX.ls...")
  error <- abs(-4125.055726045400661 - 
                est.VARX.ls(example.data)$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("       to.SS...")
  SSmodel <- l(to.SS(VARmodel),example.data)
  error <- abs(VARmodel$estimates$like[1]-SSmodel$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("     to.ARMA...")
  ARMAmodel <- l(to.ARMA(SSmodel),example.data)
  error <- abs(VARmodel$estimates$like[1]-ARMAmodel$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat(" McMillan.degree...")
  ok <- 12 == McMillan.degree(VARmodel$model, verbose=F)$distinct
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat(" Working Paper 93-4 comparisons:\n")
  if (verbose) cat("      VAR model likelihood...")
  sub.sample <- tfwindow(example.data,end=c(1981,2))
  VARmodel <- est.VARX.ar(sub.sample, re.add.means=F)
  error <- abs(-2567.3280114943772787 - VARmodel$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("      VAR model roots     ...")
  error <- abs(4.6786105186422091151 - sum(Mod(roots(VARmodel, by.poly=T))))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }


  if (verbose) cat("      SS  model likelihood...")
  SS1.model <- l(balance.Mittnik(to.SS(VARmodel), n=9),sub.sample)
  error <- abs(-2567.328011494376824 - SS1.model$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("      SS  model roots     ...")
  error <- abs(4.6786105186422082269 - sum(Mod(roots(SS1.model))))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("     ARMA model likelihood...")
  ARMA.model<- l(to.ARMA(SS1.model),sub.sample)
  error <- abs(-2567.328011494376824 - ARMA.model$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (verbose) cat("     ARMA model roots     ...")
  z <- roots(ARMA.model, fuzz=1e-4, by.poly=T)
  # the tolerance of this comparison had to be reduced because of changes
  #  from Splus 3.2 to Splus 3.3
  error <- abs(4.6786108692196117786 - sum(Mod(z)))
  ok <- 100*fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n")  }

  if (summary) 
    {cat("All example tests completed ")
     if (all.ok) cat("OK\n")
     else cat(",some FAILED! max. error = ", max.error,"\n")
    }
  invisible(all.ok)
}

example.BOC.93.4.paper.tests <- function(example.data, example.raw.data, verbose=T, fuzz.small=1e-10)
{# Use as follows:
 #   example.BOC.93.4.paper.tests(eg1.DSE.data.diff, eg1.DSE.data) 
 cat("This function reproduces some results from Bank of Canada Working Paper 93-4.\n")
   if (!example.verify.data(example.data, fuzz.small=fuzz.small))
     stop("example.data does not verify correctly. Testing stopped.")
   sub.sample <- TSdata(
      input=tfwindow(input.data(example.data),end=c(1981,2)),
      output=tfwindow(output.data(example.data),end=c(1981,2)) )
   VAR.model <- est.VARX.ar(sub.sample, re.add.means=F)
   SS1.model <- l(balance.Mittnik(to.SS(VAR.model), n=9),sub.sample)
   g1 <- diag(1,9)
   g1[1:3,] <- SS1.model$model$H
   g1 <- solve(g1)
   g2 <- diag(1,9)
   g2[3,2:3] <- c(.1,-.1)
   g2[9,9] <- -1  # this is not really necessary but seems to have
                  #   happened in the paper
   example.gap.matrix <-g1 %*% g2
   SSgap.model <- l(gmap(example.gap.matrix,SS1.model),sub.sample)
   ARMA.model<- l(to.ARMA(SS1.model),sub.sample)

# the model parameters could be displayed at this point by:
#         display(VAR.model)   
#         display(SS1.model)   etc.

   cat("Likelihood of VAR model:                          ")
   print(VAR.model$estimates$like[1], digits=16)
   cat("Likelihood of Mittnik balanced state space model: ")
   print(SS1.model$estimates$like[1], digits=16)
   cat("Likelihood of state space `gap' model:            ")
   print(SSgap.model$estimates$like[1], digits=16)
   cat("Likelihood of ARMA model:                         ")
   print(ARMA.model$estimates$like[1], digits=16)
   cat("Remark: A small change has been made in the likelihood\n")
   cat("calculation since the version of the code used for\n")
   cat("calculating the results in Bank of Canada Working Paper 93-4.\n")
   cat("The new method is more robust to degenerate densities but gives\n")
   cat("a small difference in the likelihood value. (The value reported \n")
   cat(" was -2567.32801321424.      P.Gilbert.\n")
   
   cat("Stability of VAR model:\n")
   stability(VAR.model)
   cat("Stability of Mittnik balanced state space model:\n")
   stability(SS1.model)
   cat("Stability of state space `gap' model:\n")
   stability(SSgap.model)
   cat("Stability of ARMA model:\n")
   stability(ARMA.model)

   if( dev.cur() != 1 ) 
     {tfplot(VAR.model, Title="VAR model")
      cat("Remark: These are not advertised as best estimates. There is a bias.\n")
      cat("This estimation technique may be improved by setting some of the\n") 
      cat("options and other available estimation techniques work better.\n")
      cat("The example is intended primarily for illustrating the equivalence.\n")
      cat("press return to continue>");key<-dsescan(what="");cat("\n")
      tfplot(SS1.model, Title="Mittnik balanced state space model")
      cat("press return to continue>");key<-dsescan(what="");cat("\n")
      tfplot(SSgap.model, Title="State space `gap' model")
      cat("press return to continue>");key<-dsescan(what="");cat("\n")
      tfplot(ARMA.model,  Title="ARMA model")
      cat("press return to continue>");key<-dsescan(what="");cat("\n")
      model<- l(VAR.model,example.data)  # full sample
      example.show.ytoy.cpi(model,example.raw.data)  
        title(main="Predicted and actual CPI in terms of per cent change over 12 months")
      example.show.ytoy.cpi(model,example.raw.data, start=240)    
        title(main="Predicted and actual CPI in terms of per cent change over 12 months - ex post period")
     }
   
  invisible() 
}

example.show.ytoy.cpi <- function(model, raw.data, start = 1)
{
# plot cpi in year over year % change.
# prediction is relative to previous month's actual (raw.data)
# and % change is relative to actual.
# start is the starting point for plotting.
# base is the start value of the undif, un logged series.
        i <- 3 # cpi is the third variable
	base <- raw.data$output[1, i]
	pred <- model$estimates$pred[, i]
	y <- model$data$output[, i]
	y <- cumsum(c(log(base), y))
	pred <- c(log(base), pred)	# cumsum using pred relative to actual
	pred[2:length(pred)] <- pred[2:length(pred)] + y[1:(length(pred) - 1)]
	pred <- exp(pred)
	y <- exp(y)
	pred <- 100 * ((pred[13:length(pred)] - y[1:(length(y) - 12)])/y[1:(
		length(y) - 12)])
	y <- 100 * ((y[13:length(y)] - y[1:(length(y) - 12)])/y[1:(length(y) - 
		12)])
	tfplot(tfwindow(y, start=start),tfwindow(pred, start=start)) # tsplot
             invisible()
}

example.truncate.data <- function(d.all)
{ # truncate sample  to 240 periods.
  d <- list( input= input.data(d.all)[1:240,, drop=F], 
            output=output.data(d.all)[1:240,]) 
  dimnames(d$input) <- list(NULL, "R90")
  dimnames(d$output) <- list(NULL, c("M1","GDPl2", "CPI"))
  TSdata(d)
}



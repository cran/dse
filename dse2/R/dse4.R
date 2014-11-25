#   2000/04/18 11:15:49
############################################################################

# Functions in this file are mainly for evaluating the information <<<<<<<<<<<<<
#    content of data series for predicting other series.           <<<<<<<<<<<<<


############################################################################
#
#  methods for "forecast.cov.estimators.wrt.data.subsets", "forecast.cov"
#
############################################################################


is.forecast.cov.estimators.wrt.data.subsets <- function(obj) 
{inherits(obj,"forecast.cov.estimators.wrt.data.subsets")}

print.forecast.cov.estimators.wrt.data.subsets <- function(x, 
                 digits=options()$digits)
{cat("forecast.cov.estimators.wrt.data.subsets\n")
 cat("essential data:", x$essential.data, "\n")
 cat("considering:", output.series.names(x$all.data), 
                      input.series.names(x$all.data), "\n")
 invisible(x)
}

summary.forecast.cov.estimators.wrt.data.subsets <- function(object)
  {classed(list( dseclass(obj),  #summary constructor
        horizons=object$horizons, 
        essential.data=object$essential.data,
        output.names=output.series.names(object$all.data), 
        input.names =input.series.names(object$all.data),
        nxt=NextMethod("summary")), 
    "summary.forecast.cov.estimators.wrt.data.subsets")
  }

print.summary.forecast.cov.estimators.wrt.data.subsets <- function(x,
                 digits=options()$digits)
  {cat("class: ", x[[1]], "   ")
   cat(x$horizons, " horizons\n")
   cat("essential data:", x$essential.data, "\n")
   cat("considering:",    x$output.names,  x$input.names, "\n")
   print( x$nxt )
   invisible(x)
  }


############################################################################
#
#  methods for generating test data
#
############################################################################

gen.mine.data <- function(umodel, ymodel, uinput=NULL, sampleT=100, 
	unoise=NULL, usd=1,ynoise=NULL, ysd=1, rng=NULL)
{if(!is.TSm.or.em(umodel)) TS.error.exit()
 if (is.TSestModel(umodel)) umodel <- umodel$model
 if(!is.TSm.or.em(ymodel)) TS.error.exit()
 if (is.TSestModel(ymodel)) ymodel <- ymodel$model
 if (input.dimension(ymodel) != output.dimension(umodel))
   stop("umodel output dimension must equal ymodel input dimension.")
 
 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 
 input <- output.data(simulate(umodel, input=uinput, sampleT=sampleT, 
                   noise=unoise, sd=usd))
 r <- TSdata(input  = input,
             output = output.data(simulate(ymodel, input=input,
	                             sampleT=sampleT, noise=ynoise, sd=ysd)) )
 r$umodel  <- umodel
 r$ymodel  <- ymodel
 r$uinput  <- uinput
 r$sampleT <- sampleT 
 r$unoise  <- unoise
 r$usd     <- usd
 r$ynoise  <- ynoise
 r$ysd     <- ysd
 r$rng     <- rng
 r
}

build.input.models <- function(all.data, max.lag=NULL)
{# make a list of univariate models, one for each series in input.data(data)
 #   for use by build.diagonal.model. 
 n <- input.dimension(all.data)
 multi.models <- vector("list", n)
 for (i in seq(n))
   {data <-trim.na(TSdata(output= input.data(all.data, series=i)))
    multi.models[[i]] <- TSmodel(est.VARX.ls(data, max.lag=max.lag))
   }
 multi.models
}

build.diagonal.model <- function(multi.models)
{# build one diagonal model from a list of models as returned  by 
 # build.input.models. Uses the AR part only. This can be used by gen.mine.data.
 n <- length(multi.models)
 lag <- 0
 for (i in seq(n)) lag <- max(lag, dim(multi.models[[i]]$A)[1])
 p <- 0
 for (i in seq(n))  p  <- p + dim(multi.models[[i]]$A)[3]
 model <- array(0, c(lag, p,p))
 p <- 0
 for (i in seq(n))
   {pi <- dim(multi.models[[i]]$A)[3]
    li <- dim(multi.models[[i]]$A)[1]
    model[ 1:li, (p+1):(p+pi), (p+1):(p+pi)] <-  multi.models[[i]]$A
    p <- p + pi
   }
 ARMA(A= model, B=array(diag(1, p), c(1,p,p)))
}



############################################################################
#
#  methods for stepwise mining
#
############################################################################


plot.mine.stepwise <- function(x)
  {cases <- length(x$stepwise$rss)
   o <- rev(order(x$stepwise$rss))
   vo <- dim(x$s.output.indicator)[2]
   plto <- t(matrix(1:vo, vo, cases)) * x$s.output.indicator[o,]
   if (!is.null(x$s.input.indicator))
     {vi <- dim(x$s.input.indicator)[2]
      plti <- t(matrix(-1:(-vi), vi, cases)) * x$s.input.indicator[o,]
      plt <- cbind(plti,plto)
     }
   plt[plt==0] <- NA
   matplot(0:(cases-1), plt, type="p", pch="+")
   y <- NULL
   io <-   x$io.indicator   & (1==x$lag.indicator)
   if (any(io)) y <- c(y,paste("output",  x$v.indicator[io]))
   io <-  (!x$io.indicator) & (0==x$lag.indicator)
   if (any(io)) y <- c(y,paste(" input", x$v.indicator[io]))
   cat("y axis above zero (outputs) and below zero (inputs) indicate", y, "\n")
   invisible()
  }


mine.stepwise <- function(data, essential.data=1,
      method="efroymson", f.crit=2, intercept=T,
      subtract.means=F,  standardize=F, 
      lags.in=6, lags.out=6, trend=F, plot.=T) 
{  data <- freeze(data)
   m <- ncol(input.data(data))
   p <- ncol(output.data(data))
   if(is.null(m))  m <- 0
   N <- nrow(output.data(data))
   if (standardize)
     {svd.cov <- svd(var(output.data(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, list(output=scalefac))
     }
   if (subtract.means)
    {if(m!=0)
       {input.means<-apply(input.data(data),2, mean)
        input.data(data)<-input.data(data)-t(matrix( input.means, m,N))
       }
     output.means <- apply(output.data(data),2, mean)
     output.data(data)  <- output.data(data) - t(matrix(output.means, p,N))
    }
 # The matrix Past is blocks of data:
 #  [ out-1 | out-2 | ... | out-max.lag | in | in-1 | ... | in-max.lag ]
 # so the result M has a corresponding structure.
 # If there is an input variable (m!=0) then it is shifted to give feedthrough 
 #    in one period. If lags.in=lags.out this has the result that a data point
 #    is lost at the beginning of the input series.
   if(m==0)
     {Past <- matrix(NA,N-lags.out, p*lags.out)
      io.indicator <- c(rep(T, p*lags.out))
      v.indicator  <- c(rep(1:p, lags.out)) 
      lag.indicator  <- c(t(matrix(1:lags.out, lags.out,p))) 
      for (i in 0:(lags.out-1)) 
         Past[,(1+p*i):(p*(1+i))] <-output.data(data)[(lags.out-i):(N-i-1),]
      Present <- output.data(data)[(lags.out+1):N, essential.data]
     }
   else 
     {shift <- max(lags.in+1, lags.out) # start pt. for Present
      Past <- matrix(NA,N-shift+1, p*lags.out+m*(1+lags.in))
      io.indicator <- c(rep(T, p*lags.out), rep(F, m*(1+lags.in)))
      v.indicator  <- c(rep(1:p, lags.out), rep(1:m, (1+lags.in))) 
      lag.indicator<- c(t(matrix(1:lags.out, lags.out,  p)), 
                        t(matrix(0:lags.in, (1+lags.in),m))) 
      for (i in 0:(lags.out-1)) 
        Past[,(1+p*i):(p*(1+i))] <-output.data(data)[(shift-1-i):(N-1-i),]
      for (i in 0:lags.in) 
        Past[,(p*lags.out+1+m*i):(p*lags.out+m*(1+i))] <-
                                   input.data(data) [(shift-i):(N-i),]
      Present <- output.data(data)[shift:N, essential.data]
     }
   dimnames(Past) <- list(NULL, c(
         paste(c(paste("out.v", matrix(1:p, p, lags.out), "L",sep="")),
               c(t(matrix(1:lags.out,  lags.out, p))), sep=""),
         paste(c(paste("in.v", matrix(1:m, m, lags.in), "L",sep="")),
               c(t(matrix(0:lags.in, 1+lags.in,  m))), sep="")) )
   plot. <- plot. & exists.graphics.device()
   if (plot.) par(mfcol=c(2,1))
   M <- stepwise(Past,Present, method=method,f.crit=f.crit, intercept=intercept,
                 plot=plot.)
   # Now construct an inidicator (s.indicator) of the series which are used in
   # each element of rss returned by tepwise.
   # The trick is to collapse obj$stepwise$which using obj$v.indicator so that
   #   any lags of a variable get lumped together.
   p <- v.indicator * io.indicator
   # part of the following (not the outer part) is an inner 
   # prod. with | in place of + and  & in place of *
   s.output.indicator <-  0 != (M$which %*% outer(p, 1:max(p),"==") )
   m <- v.indicator * !io.indicator
   if (max(m) !=0)
      s.input.indicator <- 0 != (M$which %*% outer(m, 1:max(m), "==") )

   M <- classed(list(stepwise=M, io.indicator=io.indicator,  #constructor
             v.indicator=v.indicator,
             lag.indicator=lag.indicator, Past=Past,
             lags.in=lags.in, lags.out=lags.out,
             s.input.indicator=s.input.indicator, 
             s.output.indicator=s.output.indicator), "mine.stepwise")
   if (plot.) plot(M)
   invisible(M)
}


############################################################################
#
#  methods for mining by splitting sample for estimation and forecast error
#
############################################################################


permute <- function(M)
  {if (is.null(M)) return(NULL)
   if (M==1) return(matrix(0:1, 2,1))
   if (M==-1) return(-matrix(0:1, 2,1))
   if (M==0) return(NULL)
   r <- permute(abs(M)-1) 
   sign(M)*rbind(cbind(r,abs(M)), cbind(r,0))
  }



mine.strip <- function(all.data, essential.data=1, 
                       estimation.sample=.5, 
                       discard.before=1, horizons=1:12,quiet=F,
                       estimation.methods=NULL,
                       step.size=NULL)
{# Calculate the predictions cov for essential.data of models estimated 
 # with estimation methods indicated by estimation.methods. 
 # estimation.methods is a list with syntax similar to programs
 #  for comparing estimation methods (eg. estimate.models), BUT ONLY 
 #  THE FIRST element (estimation method) is considered.
 # Essential.data indicates the subset of output variables to included in all
 #  models. It should be a vector of the indices. All possible combinations of
 #  input series and other output series data are considered. If omitted,
 #  essential.data is taken as the 
 #  first output series.
 # Only forecast covariances for essential data are returned.
 # discard.before is an integer indicating 1+the number of points in the
 # beginning of predictions to discard for calculating prediction covariances.
 # estimation.sample indicates the portion of the data to use for estimation.
  #If estimation.sample is an integer then it is used to indicate the number
  # of points in the sample to use for estimation. If it is a fracton it is
  # used to indicate the portion of points to use for estimation. The remainder
  # of the sample is used for evaluating predictions (ie. forecast covariance).

 # If step.size is NULL then all possible data permutations are attempted.
 #  Because S has a hard-coded limit in the number of synchronize calls this is
 #  not always possible (For loops call synchronize.) An error message:
 #    Error in synchronize(1): No room in database table
 #  If step.size is not NULL it should be a positive integer. In this case 
 #  variable permutions are divided up into
 #  steps of the given size. The result returned by the function can be used
 #  to continue from the last step:
 #      intermediate.result <- mine.strip(data, ...)
 #      intermediate.result <- mine.strip(intermediate.result)
 #      intermediate.result <- mine.strip(intermediate.result)
 #      result <- mine.strip(intermediate.result)
 #  This can be done either interactively or in a batch process, but cannot be
 #  done in a function because the database table is not cleared until the top
 #  level expression is complete.
 #  The class of an intermediate result is mine.strip.intermediate.result and
 #  the class of the final result is
 #         c("forecast.cov.estimators.wrt.data.subsets", "forecast.cov")
 #  If the final result is used in a call to mine.strip then it is just 
 #  returned, so extra calls do not cause errors and are very quick.
 #  This is useful when you are too lazy to calculate the exact number of steps.

  if (dseclass(all.data)[1] == "forecast.cov.estimators.wrt.data.subsets")
       {cat("done.\n")
        return(all.data)
       }
  if (dseclass(all.data)[1] == "mine.strip.intermediate.result")
    {r <- all.data$forecast.cov
     start <- 1+all.data$end
     estimation.sample <- all.data$estimation.sample
     discard.before <- all.data$discard.before
     quiet <- all.data$quiet
     step.size <- all.data$step.size
     variable.index <- all.data$variable.index
     m <- all.data$m
     p <- all.data$p
     multi.model <- all.data$multi.model
     essential.data <- all.data$essential.data
     estimation.methods <- all.data$estimation.methods
     horizons <- all.data$horizons
     all.data <- all.data$all.data
    }
  else
    {start <- 1
     if (estimation.sample < 1.0 )  estimation.sample <- 
           as.integer(round(estimation.sample*periods(all.data)))
     discard.before <- discard.before+estimation.sample

    #first  gen. combinations of non-essential data
     p <- output.dimension(all.data)
     m <-  input.dimension(all.data)
     M <- permute(m + p - length(essential.data) )
     # now combine essential data and permutations of non-essential data.
     if (is.null(M))
       {variable.index<-matrix(essential.data,length(essential.data),1)
        warning("essential.data seems to include all series, which does not make sense in call to mine.strip.")
       } 
     else
       variable.index<-cbind(
             t(matrix(essential.data,length(essential.data), nrow(M))), 
             array(c(0,seq(p)[-essential.data],seq(m))[1+M], dim(M)))
     r <- NULL
     if   (is.null(step.size)) step.size <- nrow(variable.index)
     else if (0 == step.size)  step.size <- nrow(variable.index)
     multi.model <- NULL
    }
 end <- min(nrow(variable.index), start+step.size-1)
 for (i in start:end )
   {data<-TSdata(output=output.data(all.data, series=variable.index[i,1:p]), 
                  input= input.data(all.data, series=variable.index[i,(p+1):(p+m)]))
    if(0==length(output.data(data))) 
      stop("The variable selection has mistakenly eliminated all output variables.")
    models <-estimate.models(data, estimation.sample=estimation.sample,
                       trend=T,quiet=quiet,
                       estimation.methods=estimation.methods)
    multi.model <- append(multi.model, list(models$multi.model[[1]]))
    rn <- forecast.cov( models$multi.model[[1]], data=data, 
                    horizons=horizons, discard.before=discard.before)
    r<- append(r, list(
           rn["forecast.cov"][[1]][[1]][,essential.data,essential.data,drop=F]))
   }
  r<-list(forecast.cov=r,all.data=all.data, essential.data=essential.data,
         variable.index=variable.index,
         estimation.methods=estimation.methods,
         multi.model=multi.model,
         horizons=horizons, 
         discard.before=discard.before)
  if (end == nrow(variable.index))
    dseclass(r) <- c("forecast.cov.estimators.wrt.data.subsets", "forecast.cov")
  else
    {r<-classed(append(r, list(estimation.sample=estimation.sample, #constructor
             quiet=quiet, step.size=step.size, end=end, m=m,p=p)),
           "mine.strip.intermediate.result" )
    }
  r
}

# z <-mine.strip(eg1.DSE.data.diff, essential.data=1, 
#      estimation.methods= list(est.VARX.ar=list(max.lag=3))) 


min.forecast.cov <- function(obj, select.series=1, verbose=T)
  {#obj is an object as returned by mine.strip
   #select the min cov for select.series only!!! at each horizon and print
   # the returned result is a vector indicating the element of forecast.cov which
   # was the min at each horizon. It is suitable as an argument to plot eg:
   #     tfplot(obj, select.cov=min.forecast.cov(obj))
   # The results of this are similar to the default results of 
   #   select.forecast.cov() cov info and information about the horizon
   #   where the model is optimal are given.

   N <- length(obj$forecast.cov)
   z <- matrix(0,length(obj$horizons),N)
   for (j in 1:N) z[,j]<-obj$forecast.cov[[j]][,select.series,select.series]
   m <- apply(z,1, min)
   r <- rep(NA,length(obj$horizons))
   for (j in 1:length(obj$horizons))
      r[j] <- (seq(N)[ z[j,]== m[j] ])[1] # only the first if more than 1 min
   if (verbose)
     {cat("              model     cov          using data\n")
      for (j in 1:length(obj$horizons))
         cat("horizon ", j,"   ",  r[j],"   ", m[j],  "   ", 
             obj$variable.index[r[j],],"\n")
     }
   invisible(r)
  }


select.forecast.cov <- function(obj, select.series=1, 
    select.cov.best=1,
    select.cov.bound=NULL,
    ranked.on.cov.bound=NULL,
    verbose=T)
  {#obj is an object as returned by mine.strip
   #select models with forecast cov for select.series meeting criteria.
   # The default select.cov.best=1 selects the best model at each horizon.
   #  select.cov.best=3 would select the best 3 models at each horizon.
   #     tfplot(select.forecast.cov(obj, select.cov.best=3))
   # If select.cov.bound is not NULL then  select.cov.best is ignored and
   #  any model which is better than the bound at all horizons is selected.
   #  select.cov.bound can be a vector of the same length as select.series,
   #  in which case corresponding elements are applied to the different series.
   #  any model which is better than the bound at all horizons is selected.
   # ranked.on.cov.bound is is used if it is not NULL and select.cov.bound is
   #  NULL. In this case select.cov.best is ignored.
   #  ranked.on.cov.bound should be a positive integer. The forecast
   #  covariances are ranked by there maximum over the horizon and the
   #  lowest number up to ranked.on.cov.bound are selected. This amounts
   #  to adjusting the covariance bound to allow for the given number of
   #  models to be selected. If select.series is a vector the results are 
   #  the best up to the given number on any series!
   # select.cov.bound can be a vector of the same length as select.series,
   #  in which case corresponding elements are applied to the different series.
   # If verbose=T then summary results are printed.
   # The returned result is a forecast.cov object like obj, but filtered
   #  to remove models which do not meet criteria.
   #     tfplot(select.forecast.cov(obj, select.cov.bound=20000))

   N <- length(obj$forecast.cov)
   r <- NULL
   if (!is.null(select.cov.bound))
     if (1 == length(select.cov.bound)) 
       select.cov.bound <- rep(select.cov.bound, length(select.series))
   for (i in 1:length(select.series)) 
     {z <- matrix(NA,length(obj$horizons),N)
      for (j in 1:N) 
         z[,j]<-obj$forecast.cov[[j]][,select.series[i],select.series[i]]
      if (!is.null(select.cov.bound))
         r <- c(r, seq(N)[apply((z <= select.cov.bound[i]),2, all)])
      else if (!is.null(ranked.on.cov.bound))
         r <- c(r, order(apply(z,2,max))[1: ranked.on.cov.bound])
      else
        {#r <- c(r, apply(z,1,sort.list)[ select.cov.best,])
         r <- c(r, apply(z,1,order)[ select.cov.best,])
        }
     }
   if (0==length(r)) stop("No forecasts meet the specified criterion.")
   r <- r[!apply(outer(r,r,"==") & 
          outer(seq(length(r)),seq(length(r)),"<"),  2,any)] #eliminate repeats
   r <- sort(r) 
   pred  <- vector("list",length(r))
   model <- vector("list",length(r))

   for (j in 1:length(r))
       {pred[[j]]           <- obj$forecast.cov[[r[j] ]]
#       model[[j]]          <- obj$multi.model[[r[j] ]]
       }
   obj$forecast.cov <- pred
#  obj$multi.model <- model
   obj$variable.index <- obj$variable.index[r,, drop=F]
   obj$selection.index <- r
   if (verbose)
     {cat("    model  using subset data series (output | input)\n")
      for (j in 1:length(obj$forecast.cov))
         cat( j,"   ", r[j],  "   ", 
             obj$variable.index[j,],"\n")
     }
   invisible(obj)
  }


exclude.forecast.cov <- function(obj, exclude.series=NULL)
  {# exlude results which depend on the indicated series from a 
   #  (forecast.cov.estimators.wrt.data.subsets forecast.cov) object.
   if (!is.null(exclude.series))
     {include<- !apply(0 != obj$variable.index[,exclude.series, drop=F], 1,any)
      obj$forecast.cov   <- obj$forecast.cov[include]
      obj$variable.index <- obj$variable.index[include,]
      obj$multi.model    <- obj$multi.model[include]
      # note all.data is not changed and variable.index still refers to it.
     }
   invisible(obj)
  }

############################################################################
#
#       procedure for testing functions     <<<<<<<<<<<<<
#
############################################################################



dse4.function.tests <- function(verbose=T, synopsis=T, 
		fuzz.small=1e-14, fuzz.large=1e-7, graphics=T)
{max.error <- 0
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

 if (synopsis & !verbose) cat("All dse4 tests ...") 
 if (verbose) cat("dse4 test 1 ... ")
  z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
  ok <- is.forecast.cov.estimators.wrt.data.subsets(z)
  all.ok <-  ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse4 test 2 ... ")
  z1 <- z$multi.model[[
       select.forecast.cov(z, select.cov.best=1, verbose=F)$selection.index[2]]]
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
  ymodel$C <- array(0, c(dim(z)[1:2], output.dimension(umodel))) 
  ymodel$C[1:(dim(z)[1]), 1:(dim(z)[2]), 1:(dim(z)[3])] <- z 
  sim.data <- gen.mine.data(umodel, ymodel,
    rng= list(kind="default",seed=c(21,46,16,12, 51, 2, 31, 8, 42, 60, 7, 3)) )
  m.step <- mine.stepwise(sim.data, method="backward")
  error <- max(abs(m.step$stepwise$rss[c(1,27)] -
               c(47.537312899054931847, 4088283.2706551752053)))
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
  invisible(all.ok)
}

dse4.graphics.tests <- function(verbose=T, synopsis=T,  pause=F)
{ if (synopsis & !verbose) cat("dse4 graphics tests ...")
  if (verbose) cat("  dse4 graphics test 1 ...")

  # If no device is active then write to postscript file 
  if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
  if(pause) dev.ask(ask=T)

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


############################################################################
#
#       end
#
############################################################################

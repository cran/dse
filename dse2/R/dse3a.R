
############################################################################

# Functions in dse3a.s and dse3b.s file are mainly for evaluating estimation techniques.

# The first group are for generating simulations (ie- generate multiple
#   stochastic simulations of a model using simulate.)
#   These are good for looking at the stochastic properties of a model,
#   but mostly these are intended
#   only for verification purposes since other functions also generate 
#   simulations and it is usually more efficient to regenerate by setting
#   the RNG and seed than it is to save data.
#   The main function in this group is monte.carlo.simulations().
#   The main object class in this group is "simulation"

# The second group are for analysing the convergence of estimators. This
#  is extended to functions of estimates (such as model roots). These
#  functions evaluate a single given estimation method with multiple
#  simulated data sets from a given "true" model.
#  The main function in this group is eval.estimation().
#  The main object classes in this group are
#    "estimation.evaluation"
#     c("roots.ee","estimation.evaluation")
#     c("TSmodel.ee","estimation.evaluation")
#     c("TSestModel.ee","estimation.evaluation")
#     c("parms.ee","estimation.evaluation")           and
#     c("roots.ee","estimation.evaluation")


# The third group applies multiple estimation techniques to a given data set.
#  This is primarily a utility for other functions 
#  The main (only) function in this group is estimate.models().
#  It returns an object of class c("estimated.models")


# The fourth group looks at the forecasts and the covariance of forecasts 
#   for multiple horizons. 
#  The simplest case horizon.forecasts() which calculates the forecast for 
#    different horizons at all periods in the sample. This is primarily
#    a utility for calculating the forecast error.
#  The next case is is estimators.horizon.forecasts() which 
#   is an extention of horizon.forecasts(). It takes specified data 
#   and estimation techniques and calculates forecasts from the estimated 
#   models.
#  The generic function forecast.cov() which considers mulitple  
#   models and calculates the cov of forecasts relative to a given data set.
#   It takes a list of models (+ trend, zero)  and calculates the cov 
#   of predictions. It uses forecast.cov.single.TSmodel()
#  The next case, forecast.cov.estimators.wrt.data() uses a list of estimation
#   methods to estimate a list of models and calculate the cov of predictions 
#   relative to one given data set.

#   The next case forecast.cov.wrt.true() takes a list of models (+ trend,
#   zero)  and calculates the cov of forecasts relative to data sets 
#   simulated with a true.model.
#   The next case, forecast.cov.estimators.wrt.data simulates data and 
#   uses a list of estimation methods to estimate a list of models, then
#   calculates the cov of predictions relative to other simulated data set.
#  The main object classes in this group are
#     c("estimators.horizon.forecasts.wrt.data") # ? "horizon.forecasts")
#     "horizon.forecasts"
#     c("multi.model.horizon.forecasts","horizon.forecasts")
#     "forecast.cov"
#     c("forecast.cov.wrt.data", "forecast.cov")
#     c("forecast.cov.wrt.true", "forecast.cov")
#     c("forecast.cov.estimators.wrt.data",  "forecast.cov")
#     c("forecast.cov.estimators.wrt.true",  "forecast.cov")


# The fifth group are some experimental estimation techniques.

############################################################################
#
#       utilities  <<<<<<<<<<
#
############################################################################


minimum.startup.lag <- function(model)UseMethod("minimum.startup.lag")
minimum.startup.lag.TSestModel <- function(model)  
    {minimum.startup.lag(model$model)}

minimum.startup.lag.ARMA <- function(model)
  {lag <- dim(model$A)[1] 
   if (!is.null(model$C)) lag <- max(lag, dim(model$C)[1])
   lag
  }

minimum.startup.lag.SS <- function(model)  { 1+dim(model$F)[2] }

start.shift <- function(model,data, y0=NULL)
 {# there is some redundancy between this and  minimum.startup.lag which 
  #   should be cleaned up.
  # This function is used to determine the number of lags (and leads) needed
  #   for a model, and whether the data can be padded with zeros or the start
  #   (and end) have to be shifted within the data. Shifting is indicated if the
  #   model has an element $no.zeros (which would be specified if, for example,
  #   the model takes logrithms of data) or the data come from an external data base.
  # As of Nov. 1995 it is used by l.troll, simulate.troll, monte.carlo.troll and
  #  by feather.forecasts and forecast.cov to determine lags and whether the
  #  starting point is shifted or zeros prepended to the data.
  # if (is.TSPADIdata(data)) shift <- T requires library dsepadi
  if (inherits(data, "TSPADIdata"))  shift <- T
  else if (!is.null(model$no.zeros)) shift <-model$no.zeros
  else shift <- F

  if(!is.null(model$order))
     lags<-max(model$order$a,model$order$b, model$order$c)
  else if(!is.null(y0)) lags <- dim(y0)[1]
  else lags <-20 

  if (is.null(model$order$a.leads) & is.null(model$order$b.leads) & 
      is.null(model$order$c.leads) )      terminal.periods <- 0
  else
      terminal.periods <- max(model$order$a.leads, 
                              model$order$b.leads, model$order$c.leads)

  list(shift=shift, lags=lags, terminal.periods=terminal.periods)
 }


############################################################################
#
#       methods for monte.carlo.simulations  <<<<<<<<<<
#
############################################################################

generate.model.SS <- function(m,n,p, stable=F)
 {#randomly generate an innov state space model. Discard models with largest root
  # greater than 1 (if stable=F) or equal to or greater than 1 if stable=T.
  repeat
    {FF <- matrix(runif(n^2, min=-1,max=1),n,n)
     if (m!=0) G <- matrix(runif(n*m, min=-1,max=1),n,m)
     else G <- NULL
     H <- matrix(runif(n*p, min=-1,max=1),p,n)
     K <- matrix(runif(n*p, min=-1,max=1),n,p)
     model <- SS(F=FF, G=G,H=H,K=K)
     if (stable) {if (max(Mod(roots(model))) <  1.0) break()}
     else        {if (max(Mod(roots(model))) <= 1.0) break()}
    }
  model
 }


monte.carlo.simulations <- function( model, simulation.args=NULL, 
           replications=100, rng=NULL, ...)
{#Produces multiple simulations.
	UseMethod("monte.carlo.simulations")
}

monte.carlo.simulations.TSestModel <- function(model, simulation.args=NULL, 
           replications=100, rng=NULL, ...)
  {if (is.null(simulation.args$sd) & is.null(simulation.args$SIGMA)) 
     simulation.args$SIGMA <- model$estimates$cov
   if (is.null(simulation.args$input)) simulation.args$input <- input.data(model)
   monte.carlo.simulations(TSmodel(model), simulation.args=simulation.args, 
           replications=replications, rng=rng, ...)
  }

monte.carlo.simulations.estimation.evaluation <- function(model,...)
       {monte.carlo.simulations(TSmodel(model), rng=get.RNG(model), ...)}
monte.carlo.simulations.monte.carlo.simulation <- function(model,...)
       {monte.carlo.simulations(TSmodel(model), rng=get.RNG(model), ...)}


monte.carlo.simulations.TSmodel <- function( model, simulation.args=NULL,
          replications=100, rng=NULL, Spawn=.SPAWN, quiet=F)
{ 
 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 
 arglist <- append(list(model), simulation.args)
 if (Spawn)
  {if (!quiet)cat("Spawning processes to calculate ", replications, " replications.\n")
   assign("sim.forloop.n", replications, where = 1)
   assign("sim.forloop.result", list(NULL), where = 1)
 #  assign("sim.forloop.model", model, where = 1)
   assign("sim.forloop.arglist", arglist, where = 1)
   on.exit(remove(c("sim.forloop.i", "sim.forloop.n", "sim.forloop.result",
       "sim.forloop.arglist"),where = 1))
   For(sim.forloop.i = 1:sim.forloop.n, sim.forloop.result[[sim.forloop.i]] <- 
       do.call("simulate",  sim.forloop.arglist),
       first=options(warn=-1), sync = T)
   result <- array(NA, c(dim(sim.forloop.result[[1]]$output),replications))
   tfr <- tframe(sim.forloop.result[[1]]$output)
   for (i in 1:replications) result[,,i] <- sim.forloop.result[[i]]$output
  }
 else
  {#r <- simulate(model, list(...))$output
   r <- do.call("simulate", arglist)$output
   result <- array(NA, c(dim(r),replications))
   tfr <- tframe(r)
   result[,,1] <- r
   if (1 < replications)
     for (i in 2:replications) 
        result[,,i] <- do.call("simulate", arglist)$output
  }
series.names(result) <- output.series.names(model)
result <- tframed(result, tfr)  # my more general multidimensional ts
invisible(classed( # constructor monte.carlo.simulations
         list(simulations=result,
              model=model, rng=rng, version=version, 
              simulation.args=simulation.args,
              description = "data generated by monte.carlo.simulation.default"),
   c("monte.carlo.simulations") ))
}

is.monte.carlo.simulation <- function(obj) 
   {inherits(obj,"monte.carlo.simulations")}

print.monte.carlo.simulations <- function(x, digits=4)
{cat("Simulation with RNG ", x$rng, " from model:\n")
 print(x$model)
 invisible(x)
}

output.series.names.monte.carlo.simulations <- function(obj)
   {dimnames(obj$simulations)[[2]]}

input.series.names.monte.carlo.simulations <- function(obj)
  {input.series.names(obj$simulation.args$data)}

test.equal.monte.carlo.simulations <- function(d1,d2, fuzz=1e-16)
 {if (length(d1$result) != length(d2$result)) r <-F
  else  r <- all(fuzz > abs(d1$simulations - d2$simulations))
  r
 }


summary.monte.carlo.simulations <- function(object,
        select.series=NULL, periods=1:3)
 {stats <- NULL
  if (!is.null(select.series))
    {if (sim.dim[3] <20) 
        warning("SD calculation is not very good with so few simulations.")
     names <- output.series.names(object)
     if(is.null(names)) names <- output.series.names(object$model)
     if (!is.numeric(select.series)) select.series <-match(select.series, names)
     names <- names[select.series]
     mn<-apply(object$simulations[periods,select.series,,drop=F],c(1,2),mean)
     sd<-apply(object$simulations[periods,select.series,,drop=F],c(1,2),var)^0.5
     stats <- rbind(mn,sd) 
     dimnames(stats)<- list(c(paste("mean period", periods), 
                          paste("S.D. period",periods)), names)
    }
  classed(list(  # constructor summary.monte.carlo.simulations
     description=object$description,
     sampleT=dim(object$simulations)[1],
     p=      dim(object$simulations)[2], 
     simulations=dim(object$simulations)[3], 
     summary.stats=stats,
     rng=get.RNG(object)), 
  "summary.monte.carlo.simulations")
}


print.summary.monte.carlo.simulations <- function(x, digits=options()$digits)
 {cat("Object class monte.carlo.simulations\n")
  cat(x$description, "\n")
  cat("periods=",x$sampleT, "variables=", x$p,"simulations=",x$simulations,"\n")
  cat("rng= ", x$rng, "\n")
  if (!is.null(x$summary.stats))   print(x$summary.stats, digits=digits)
  invisible(x)
 }



tfplot.monte.carlo.simulations <- function(x, start.=NULL, end.=NULL,
    select.series=seq((dim(x$simulations)[2])), 
    select.simulations=seq(dim(x$simulations)[3]),
    graphs.per.page=5)
  {names <- output.series.names(x)
   if(is.null(names)) names <- output.series.names(x$model)
   if (is.null(start.)) start. <- start(x$simulations)
   if (is.null(end.))   end.   <- end(x$simulations)
   tf.p <- tframe(x$simulations)
   Ngraph <- min(length(select.series), graphs.per.page)
   old.par <-par(mfcol = c(Ngraph, 1), mar= c(5.1,6.1,4.1,2.1))
   on.exit(par(old.par))
   #zz<- matrix(NA, dim(sim)[1], length(x$simulations))
   if (!is.numeric(select.series)) select.series <- match(select.series, names)
   for(i in select.series) 
        {zz <- (x$simulations)[,i,select.simulations]
         tframe(zz) <- tf.p
         tfplot(tfwindow(zz,start=start.,end=end., warn=F), ylab=names[i]) #tsplot
         if(i == select.series[1])  title(main = "Monte Carlo Simulations")
        }
   invisible()
}



distribution.monte.carlo.simulations <- function(obj,
     select.series=seq(dim(obj$simulations)[2]),
     x.sections=T, periods=1:3, graphs.per.page=5)
  {
if (dim(obj$simulations)[3] <20) 
     warning("This is not very good with so few simulations.")
   names <- output.series.names(obj)
   if(is.null(names)) names <- output.series.names(obj$model)
   if (!is.numeric(select.series)) select.series <- match(select.series, names)
   names <- names[select.series]
   Ngraph <- min(length(select.series), graphs.per.page)
   if (x.sections)
       {data <- obj$simulations[periods, select.series,,drop=F]
        old.par <-par(mfrow =c(Ngraph, length(periods)), mar=c(5.1,6.1,4.1,2.1))
       }
   else 
       {old.par <-par(mfrow =c(Ngraph, 1), mar=c(5.1,6.1,4.1,2.1))
        mn <- apply(obj$simulations[, select.series,,drop=F], c(1,2), mean)
        sd <- apply(obj$simulations[, select.series,,drop=F], c(1,2), var) ^ 0.5
        plt <- array(c(mn, mn+sd, mn-sd, mn+2*sd, mn-2*sd), c(dim(mn),5))
        tf.p <- tframe(obj$simulations)
      }
   on.exit(par(old.par))
   for (i in 1:length(select.series)) 
    {if (x.sections)
        {for (j in 1:length(periods)) 
          {if (exists("ksmooth")) plot(ksmooth(data[j,i,], # -mean[j,i],
                              bandwidth=var(data[j,i,])^0.5, kernel="parzen"),
                        type="l",xlab=paste("Period",periods[j]),ylab=names[i])
           else if (exists("density")) plot(density(data[j,i,]), # -mean[j,i]
                 type="l",xlab=paste("Period",periods[j]),ylab=names[i])
           else
        stop("Neither ksmooth nor density are available to calculate the tfplot.")
           if ((i == 1) & (j ==length(periods)%/%2))
              title(main = "kernel estimate of distributions")
           }
        }
     else
        {pl <-plt[,i,]
         tframe(pl) <- tf.p
         tfplot(pl, type="l", lty= c(1,3,3,2,2), ylab=names[i]) #tsplot
         if (i == 1) title(main = "Simulation mean, 1 & 2 S.D. estimates")
        }
    }
  invisible()
  }

############################################################################
#
#       methods for estimation.evaluation.  <<<<<<<<<<
#
############################################################################

#e.bb.ar.100 <- eval.estimation( mod2, replications=100, 
#               estimation.args=list(estimation="est.VARX.ar", verbose=F))

#e.bb.ls.over <- eval.estimation( simple.mod, replications=100, 
#   estimation.args=list(estimation="est.VARX.ls", max.lag=6, verbose=F), 
#   criterion="parms")


eval.estimation <- function( model, replications=100, rng=NULL, quiet=F, 
                       simulation.args=NULL,
                       estimation=NULL, estimation.args=NULL, 
                       criterion ="parms", criterion.args =NULL, spawn=.SPAWN)
{# estimation.args and criterion.args should be NULL if no args are needed.
 # If model is an object of class "estimation.evaluation" or "simulation"
 # then the model and the seed!!! are extracted so the evaluation will be
 # based on the same generated sample.
 # criterion can be in { "parms", "roots", TSmodel", "TSestModel"}
 # With the default (parms) and $model the other criteria can be reconstructed
 #  when the estimation method gets the correct form for the model. ( This
 #  is not usually the case with the default method "black.box".)
 # This is done by via the generic functions roots, TSmodel and TSestModel.
 # If criterion = "roots" then criterion.args= list(verbose=F) is advised.
 # example simulation.args=list(sampleT=100, sd=1.5)

 if(is.null(estimation)) stop("estimation method must be specified.")
 if (is.estimation.evaluation(model) | is.monte.carlo.simulation(model))
   {rng  <- get.RNG(model)
    model<- TSmodel(model)
   }
  
 truth <- do.call(criterion, append(list(model), criterion.args))

 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 
 if (!spawn)
   {if(!quiet) cat("Calculating ", replications, " estimates.\n")
    result <- vector("list",replications)
    for (i in 1:replications)
       {data <- do.call("simulate", append(list(model), simulation.args))
        m   <-  do.call(estimation, append(list(data),  estimation.args))
        result[[i]]<-do.call(criterion, append(list(m), criterion.args))
       }
   }
 else
    {if(!quiet)
        cat("Spawning processes to calculate ", replications, " estimates.\n")
     est.forloop <- function(estimation, estimation.args, model, 
                             simulation.args, criterion, criterion.args)
       {data <- do.call("simulate", append(list(model), simulation.args))
        m   <-  do.call(estimation, append(list(data),  estimation.args))
        do.call(criterion, append(list(m), criterion.args))
       }
     assign("est.forloop", est.forloop, where = 1)
     assign("est.forloop.n", replications, where = 1)
     assign("est.forloop.result", list(NULL), where = 1)
     assign("est.forloop.estimation", estimation, where = 1)
     assign("est.forloop.model", model, where = 1)
     assign("est.forloop.simulation.args", simulation.args, where = 1)
     assign("est.forloop.criterion", criterion, where = 1)
     assign("est.forloop.estimation.args", estimation.args, where = 1)
     assign("est.forloop.criterion.args", criterion.args, where = 1)

     on.exit(remove(c("est.forloop", "est.forloop.i", "est.forloop.n",
         "est.forloop.result",  "est.forloop.estimation","est.forloop.model",
         "est.forloop.simulation.args", "est.forloop.criterion", 
         "est.forloop.estimation.args","est.forloop.criterion.args"),where = 1))

     For(est.forloop.i = 1:est.forloop.n, est.forloop.result[[est.forloop.i ]]<-
       est.forloop(est.forloop.estimation, est.forloop.estimation.args, est.forloop.model, 
       est.forloop.simulation.args, est.forloop.criterion, est.forloop.criterion.args),
       first=options(warn=-1), sync = T)
     result<-est.forloop.result
    }
invisible(classed( # constructor estimation.evaluation (eval.estimations)
      list(result=result,truth=truth,model=model,
           rng=rng, version=version,
           estimation=estimation, estimation.args=estimation.args,
            criterion=criterion,   criterion.args=criterion.args, 
            simulation.args=simulation.args),
     c(paste(criterion,".ee",sep=""), "estimation.evaluation")))
}

eval.estimation.set <- function(estimation, model, 
         replications=10, sampleT=100,
         eval.proc.args=list(NULL), 
         criterion = "roots", p.set=NULL)
{ # evaluate for a set of models.   NOT WORKING
 r <- NULL
   for (p in p.set)
     {m<-model
      m$parms <- p
      m <- set.arrays(m)
      r[[p]] <- estimation.test(estimation, m, replications=10, sampleT=100,eval.proc.args=list(NULL), criterion = "roots")
     }
r
}

is.estimation.evaluation <- function(obj){inherits(obj,"estimation.evaluation")}

test.equal.estimation.evaluation <- function(obj1,obj2)
 {all(as.character(obj1) == as.character(obj2))}

print.estimation.evaluation <- function(x, digits=4)
{cat("Estimation evaluation with model:\n")
 print(x$model, digits=digits)
 cat("Evaluation criterion: ",x$criterion, "\n")
 invisible(x)
}

summary.estimation.evaluation <- function(object)
 {classed(list( # constructor summary.estimation.evaluation
     class=dseclass(object),
     estimation=object$estimation,
     estimation.args= if(!is.list((object$estimation.args)[[1]]))
                           object$estimation.args    else    NULL,
     criterion=object$criterion,
     labels=object$labels,
     criterion.args=object$criterion.args,
     replications=length(object$result), 
     true.model=object$model,
     rng=get.RNG(object)), 
  "summary.estimation.evaluation")
 }


print.summary.estimation.evaluation <- function(x, digits=options()$digits)
{ cat("Object of class: ", x$class, "\n")
  cat("Evaluation of `",x$estimation,"'")
  if(!is.list((x$estimation.args)[[1]]))
    {cat( " estimation with argument ") 
     cat(labels(x$estimation.args),"= `",x$estimation.args,"'")
    }
  cat("\n")
  cat("using criterion `", x$criterion, "' with argument ")
  cat(x$labels," = `", x$criterion.args, "'\n")
  cat(x$replications, " replications, RNG = ", x$rng, "\n")
  cat("true model:\n")
  print(x$model)
  invisible(x)
 }



distribution <- function(obj, ...)UseMethod("distribution")

distribution.TSdata <- function(obj, bandwidth=0.2, series=NULL)
  {if (0 !=  input.dimension(obj)) distribution( input.data(obj))
   if (0 != output.dimension(obj)) distribution(output.data(obj))
  }

distribution.default <- function(obj, bandwidth=0.2, series=NULL)
  {# obj should be a ts matrix (perhaps this should be a tf method).
   # If series is NULL then all series are ploted.
   # note that this graphic can be fairly misleading:
   #    distribution(runif(1000))  should be uniform
   names <- series.names(obj)
   if (!is.matrix(obj) ) obj <- matrix(obj, length(obj), 1)
   if(!is.null(series)) 
     {obj   <-   obj[,series, drop=F]
      names <- names[series]
     }
   par(mfcol=c(ncol(obj),1))
   for ( i in 1:ncol(obj))
      {if      (exists("ksmooth")) rd <- ksmooth(obj[,i], bandwidth=bandwidth) 
       else if (exists("density")) rd <- density(obj[,i], bw= bandwidth)
       else     stop("Neither ksmooth nor density are available.")
       plot(rd, type="l", ylab="density", ylim=c(0, max(rd$y)), xlab=names[i] )
      }
   invisible()
  }


distribution.estimation.evaluation <- function(obj, ...)
 {distribution(parms(obj), ...)}



############################################################################
#
#       methods for roots.ee  (estimation.evaluation)  <<<<<<<<<<
#
############################################################################

summary.roots.ee <- function(object, verbose=T)
{ nxt <- if (verbose) NextMethod("summary") else NULL
  if (! verbose) conv <- NULL
  else 
    {if (!is.null(object$result.conv)) conv <- NULL
     else                              conv <- object$result.conv
    }
  N <- length(object$result)
  p <- 0
  for (i in 1:N) p <- max(p, length((object$result)[[i]]))
  r <- matrix(NA, N, p)
  for (i in 1:N) r[i,1:length((object$result)[[i]])] <- (object$result)[[i]]
  m <- apply(r,2,sum)/N
  cov <- r- t(matrix(object$truth, p, N))
  cov <- (t(Conj(cov)) %*% cov)/(N-1)
  ecov <- r- t(matrix(m, p, N))
  ecov <- (t(Conj(ecov)) %*% ecov)/(N-1)
  classed(list( # constructor summary.summary.roots.ee
     nxt=nxt,
     conv=conv,
     true.criterion=object$truth,
     mean=m,
     cov=cov,
     ecov=ecov),
  "summary.summary.roots.ee")
 }

print.summary.roots.ee <- function(x, digits=options()$digits)
 {if (!is.null(x$nxt)) print(x$nxt)
  if (!is.null(x$conv))
        {if(all(x$conv)) cat("All estimates converged.\n")
         else cat(sum(!x$conv)," estimates did not converge!\n")
        }
  cat("\nTrue model criterion mean: ",x$true.criterion,"\n")
  cat("Sampling estimate of mean: ",x$mean,"\n")
  cat("Estimate of sampling covariance [e*Conj(t(e))] using true model:\n")
  print(x$cov)
  cat("\nEstimate of sampling covariance (using sample mean and not the true model):\n")
  print(x$ecov)
  invisible(x)
 }


tfplot.roots.ee <- function(x, ...){UseMethod("plot.roots.ee")}

plot.roots.ee <- function(x, complex.plane=T, cum=T, norm=F, bounds=T, transform=NULL, invert=F, Sort=T)
{# If complex.plane is T then all results are plotted on a complex plane and 
 #   the arguements cum and Sort do not apply. If complex.plane is F 
 #   then a sequential plot of the real and imaginary parts is produced.
 # If cum is true the cummulative average is plotted.
 # If mod is true the modulus is used, otherwise real and imaginary are separated.
 # if invert is true the reciprical is used (before cummulating).
 # if Sort is true then sort is applied (before cum but after mod) by the Re part of the root.
 #   Some grouping is usually necessary since roots are not in an obvious order
 #   but sorting by the real part of the roots could be improved upon.
   N<-length(x$result)
   n <- 0
   for (i in 1:N) n <- max(n, length((x$result)[[i]]))
   r <- matrix(0,N, n) 
   for (i in 1:N) r[i,1:length((x$result)[[i]])] <- (x$result)[[i]]
   true.lines <- c(x$truth, rep(0,n-length(x$truth)))
   if (invert)
         {true.lines <- 1/true.lines
          r <- 1/r
         }
   if(!is.null(transform)) 
         {r <- do.call(transform,list(r))
          true.lines <-do.call(transform,list(true.lines))
         }
   if (complex.plane)
    {plot.roots(x$truth, pch="o")
     for (i in 1:N) add.plot.roots(r[i,], pch="*") 
     add.plot.roots(0, pch="+") # add.plot.roots(0+0i, pch="+")
    }
  else
     {if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      if (cum) r <- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      else     r <- r 
      if(is.complex(r))
         { r <- cbind(Re(r), Im(r))
          true.lines <-c(Re(true.lines),Im(true.lines))
         }
      r[is.infinite(r)] <- 0
      true.lines <-t(matrix(true.lines, length(true.lines),N))
      matplot(x=seq(nrow(r)), y=cbind(0,true.lines, r), type="l",
             lty=c(1,rep(3,dim(true.lines)[2]), rep(2,dim(r)[2])) )
     }
  invisible(r)
}

roots.roots.ee <- function(obj, ...)   {obj}

distribution.roots.ee <- function(obj, mod=T, invert=F, Sort=F, bandwidth=0.2, select=NULL)
{# if mod is true the modulus is used, otherwise real and imaginary are separated.
 # if invert is true the reciprical is used.
 # if Sort is true then sort is applied (before cum). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
 # If select is not NULL then only the indicated roots are plotted. 
 #     ie - select=c(1,2)  will plot only the two largest roots
      N<-length(obj$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,] <- c((obj$result)[[i]], 
                                rep(0,n-length((obj$result)[[i]])))
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (invert)
         {true.lines <- 1/true.lines
          r <- 1/r
         }
      if(mod) 
         {r <- Mod(r) 
          true.lines <-Mod(true.lines)
          xlab <-"Mod root "
         }
      else
         { r <- cbind(Re(r), Im(r))
          true.lines <-c(Re(true.lines),Im(true.lines))
          xlab <-"Real part root "
         }
      r[is.infinite(r)] <- 0
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      if(!is.null(select)) r <- r[,select, drop=F]
      par(mfcol=c(dim(r)[2],1))
      for ( i in 1:dim(r)[2])
         {if      (exists("ksmooth")) rd <- ksmooth(r[,i], bandwidth=bandwidth) 
          else if (exists("density")) rd <- density(r[,i], bw= bandwidth)
          else
        stop("Neither ksmooth nor density are available to calculate the tfplot.")
          if (i > n) xlab <-"Imaginary part root "
          plot(rd, type="l", ylab="density", ylim=c(0, max(rd$y)),
               xlab=paste(xlab, n-(-i%%n)) )
          lines(rep(true.lines[i],2),c(1,0))
         }
      invisible()
}


############################################################################
#
#       methods for parms.ee (estimation.evaluation)  <<<<<<<<<<
#
############################################################################

summary.parms.ee <- function(object, verbose=T)
  {classed(summary.roots.ee(object, verbose=verbose), "summary.parms.ee")} # constructor


print.summary.parms.ee <- function(x, digits=options()$digits)
{UseMethod("print.summary.roots.ee")}


tfplot.parms.ee <- function(x, cum=T, norm=F, bounds=T, invert=F, Sort=F)
{# if cum is true the cummulative average is plotted.
 # if norm is true the norm is used, each parameter is plotted.
 # if invert is true the reciprical is used (before cummulating).
 # if Sort is true then sort is applied (before cum). This is not usually
 #   recommended but of interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
      N<-length(x$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((x$result)[[i]]))
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,1:length((x$result)[[i]])] <- (x$result)[[i]]
      true.lines <- c(x$truth, rep(0,n-length(x$truth)))
      if (invert)
         {true.lines <- 1/true.lines
          r <- 1/r
         }
      if(norm) 
         {r <- matrix((apply(r^2,1,sum))^.5, N,1)
          true.lines <-sum(true.lines^2)^.5
         }
      r[is.infinite(r)] <- 0
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      true.lines <-t(matrix(true.lines, length(true.lines),N))
      Om <-NULL
      if (bounds)
        {z  <- r-true.lines
         Om <- t(z) %*% z/(nrow(z)-1)
         Om <- diag(Om)^.5
         Om <- t(matrix(Om, length(Om), N))
         Om <- Om/matrix((1:N)^.5 , N, ncol(Om))
         Om <- cbind(true.lines+Om, true.lines-Om)
        }
      if (cum) r<- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      matplot(x=matrix(seq(nrow(r)),nrow(r),1), y=cbind(0,true.lines,r, Om), 
              type="l", lty=c(1,rep(3,dim(true.lines)[2]), rep(4,dim(r)[2]), 
                     rep(2,2*dim(r)[2]))) 
      invisible(r)
}

roots.parms.ee <- function(obj, criterion.args=NULL)
{# extract roots criterion 
  model <- obj$model
  truth <-do.call("roots", append(list(model), criterion.args))
  r <- NULL
  for (m in obj$result)
    {model$parms <- m
     model <- set.arrays(model)
     r <- append(r, 
           list(do.call("roots", append(list(model), criterion.args))))
    }
  ok <- T
  for (m in obj$result)
     ok <- ok & (length(model$parms) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("roots.ee","estimation.evaluation")))# constructor
}

distribution.parms.ee <- function(obj,  Sort=F, bandwidth=0.2)
{# if Sort is true then sort is applied (before cum). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
      N<-length(obj$result)
      n <- length(obj$truth)
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,1:length((obj$result)[[i]])] <- (obj$result)[[i]]
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      xlab <-"parameter "
      par(mfcol=c(dim(r)[2],1))
      for ( i in 1:dim(r)[2])
         {if (is.Splus()) rd <- ksmooth(r[,i], bandwidth=bandwidth)
          if (is.R())     rd <- density(r[,i], bw=bandwidth)
          plot(rd, type="l", ylim=c(0, max(rd$y)),
               ylab="density",  xlab=paste(xlab, i) )
          lines(rep(true.lines[i],2),c(1,0))
         }
      invisible()
}

TSmodel.parms.ee <- function(obj)
{# rebuild model from parms
  model <- obj$model
  truth <-TSmodel(model)
  r <- NULL
  for (m in obj$result)
    {model$parms <- m
     model <- set.arrays(model)
     r <- append(r, list(model))
    }
  ok <- T
  for (m in obj$result)
     ok <- ok & (length(model$parms) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSmodel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSmodel.ee","estimation.evaluation")))
}

TSestModel.parms.ee <- function(obj)
{# rebuild ... 
  model <- obj$model
  truth <-l(TSmodel(model), data)   # need to regenerate data
  r <- NULL
  for (m in obj$result)
    {model$parms <- m
     model <- l( set.arrays(model), data)
     r <- append(r, list(model))
    }
  ok <- T
  for (m in obj$result)
     ok <- ok & (length(model$parms) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSestModel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSestModel.ee","estimation.evaluation")))
}

############################################################################
#
#       methods for TSmodel.ee (estimation.evaluation)  <<<<<<<<<<
#
############################################################################

summary.TSmodel.ee <- function(object)
 {if (is.null((object$result)[[1]]$converged)) conv <- NULL
  else
    {conv <- rep(NA,length(object$result))
     for (i in 1:length(conv)) conv[i] <- (object$result)[[i]]$converged
    }
 # summary(parms(object))
 #summary(roots(object))  these are slow

  classed(list( # constructor summary.TSmodel.ee
     class=dseclass(object),
     conv=conv,
     default=summary.default(object)),
  "summary.TSmodel.ee")
 }

print.summary.TSmodel.ee <- function(x, digits=options()$digits)
{ cat("Object of class: ",x$class, "\n")
  if (!is.null(x$conv))
        {if(all(x$conv)) cat("All estimates converged.\n")
         else cat(sum(!x$conv)," estimates did not converge!\n")
        }
  print(x$default)
  invisible(x)
}


parms.TSmodel.ee <- function(obj, criterion.args=NULL)
{# extract parameters from models in the list and 
 #   return a list of class parms.ee estimation.evaluation 
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodel.ee can provide this argument.
  truth <-parms(obj$truth)
  r <- NULL
  for (m in obj$result) 
     r <- append(r,list(parms(m)))
  if (! is.null((obj$result)[[1]]$converged))
    {rc <- rep(NA,length(obj$result))
     for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
     obj$result.conv<-rc
    }
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"parms"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("parms.ee","estimation.evaluation")))
}

roots.TSmodel.ee <- function(obj, criterion.args=list( randomize=T))
{# extract roots criterion 
  truth <-do.call("roots", append(list(obj$truth), criterion.args))
  r <- NULL
  for (m in obj$result)
     r <- append(r, 
           list(do.call("roots", append(list(m), criterion.args))))
  if (! is.null((obj$result)[[1]]$converged))
    {rc <- rep(NA,length(obj$result))
     for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
     obj$result.conv<-rc
    }
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("roots.ee","estimation.evaluation")))
}


tfplot.TSmodel.ee <- function(x, graph.args=NULL,
                       criterion ="parms", criterion.args=NULL)
{# extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSmodel.ee", sep=""), 
               append(list(x), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
}

############################################################################
#
#       methods for TSestModel.ee (estimation.evaluation)   <<<<<<<<<<
#
############################################################################

summary.TSestModel.ee <- function(object)
  {classed(summary.TSmodel.ee(object), "summary.TSestModel.ee") }  # constructor 

print.summary.TSestModel.ee <- function(x, digits=options()$digits)
  { UseMethod("print.summary.TSmodel.ee")}


parms.TSestModel.ee <- function(obj, criterion.args=NULL)
{# extract parameters from models in the list and convergence info.
 #   return a list of class parms.ee estimation.evaluation
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodel.ee can provide this argument.
  truth <-parms(obj$truth)
  r <- NULL
  for (m in obj$result) r <- append(r,list(parms(m)))
  rc <- rep(NA,length(obj$result))
  for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
  obj$result<-r
  obj$result.conv<-rc
  obj$truth <-truth
  obj$criterion<-"parms"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("parms.ee","estimation.evaluation")))
}

roots.TSestModel.ee <- function(obj, criterion.args=NULL)
{# extract roots criterion 
  truth <-do.call("roots", append(list(obj$truth), criterion.args))
  r <- NULL
  for (m in obj$result)
     r <- append(r, 
             list(do.call("roots", append(list(m), criterion.args))))
  rc <- rep(NA,length(obj$result))
  for (i in 1:length(rc)) rc[i] <- (obj$result)[[i]]$converged
  obj$result<-r
  obj$result.conv<-rc
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("roots.ee","estimation.evaluation")))
}


tfplot.TSestModel.ee <- function(obj, graph.args=NULL,
                       criterion ="parms", criterion.args=NULL)
{# extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSestModel.ee", sep=""), 
               append(list(obj), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
}

############################################################################
#
#       function for generating estimated.models (and methods).   <<<<<<<<<<
#
############################################################################

estimate.models <- function(data, estimation.sample=NULL, trend=F,quiet=F,
                       estimation.methods=NULL)
{# Estimate models from data with methods indicated by estimation.methods. 

  if (!is.null(estimation.sample))
    {# is.integer in the next line does not work 
     if (0 != (estimation.sample %%1))
        stop("estimation.sample must be an integer.")
     if (estimation.sample <= 0)
        stop("estimation.sample must be a positive integer.")
     if (nrow(output.data(data)) < estimation.sample)
        stop("estimation.sample cannot be greater than the sample size.")
     output.data(data) <- output.data(data)[1:estimation.sample,, drop=F]
     if (0 != (input.dimension(data)))
        input.data(data) <- input.data(data)[1:estimation.sample,, drop=F]
    }
   r <-list(estimation.methods=estimation.methods)
  if (trend) r$trend.coef <- lsfit(1:periods(data), output.data(data))$coef
  if (!is.null(estimation.methods))
    {r$multi.model <- vector("list", length(estimation.methods))
     for (j in 1:length(estimation.methods))
       {est <-  do.call(names(estimation.methods)[j], 
                  append(list(data),  estimation.methods[[j]]))
        if (!is.null(est$converged) )
              est$model$converged <- est$converged
        # $ causes problems here
        r$multi.model[[j]] <- est$model
       }
    }
  classed(r, "estimated.models")
}

is.estimated.models <- function(obj){ inherits(obj,"estimated.models") }

test.equal.estimated.models <- function(obj1,obj2)
 {all(as.character(obj1) == as.character(obj2))}

print.estimated.models <- function(x, digits=4)
 {cat("Estimated models:\n")
  if (!is.null(x$trend.coef)) cat("trend coef: ", x$trend.coef, "\n")
  if (!is.null(x$multi.model))
    {for (j in 1:length(x$multi.model))
       {cat("model ", j, "\n")
        print((x$multi.model)[[j]])
        cat("\n")
       }
    }
 invisible(x)
}

summary.estimated.models <- function(object)
 {if (!is.null(object$trend.coef)) cat("trend coef: ", object$trend.coef, "\n")
  if (!is.null(object$multi.model))
    {estimation.names   <- vector(NA, length(object$multi.model))
     estimation.methods <- vector(NA, length(object$multi.model))
     for (j in 1:length(object$multi.model))
        {estimation.names[j] <- name(object$estimation.methods)[j]
         estimation.methods[j] <-    object$estimation.methods[[j]]
        }
    }
  classed(list(  # constructor summary.estimated.models
     class=dseclass(object),
     trend.coef=object$trend.coef,
     estimation.names=estimation.names,
     estimation.methods=estimation.methods,
     conv=conv),
  "summary.estimated.models")
 }

print.summary.estimated.models <- function(x, digits=options()$digits)
{ cat("Object of class: ",x$class, "\n")
  cat("Estimated models:\n")
  if (!is.null(x$trend.coef)) cat("trend coef: ", x$trend.coef, "\n")
  if (!is.null(object$multi.model))
    {for (j in 1:length(x$estimation.methods))
        cat(x$estimation.names[j],  x$estimation.methods[j], "\n")
    }
  invisible(x)
}



roots.estimated.models <- function(obj, digits=4, mod=F)
 {cat("Estimated models:\n")
  if (!is.null(obj$trend.coef)) cat("trend coef: ", obj$trend.coef, "\n")
  if (!is.null(obj$multi.model))
    {r <- vector("list",length(obj$multi.model)) 
     for (j in 1:length(obj$multi.model))
       {cat("model ", j, "\n")
        r[[j]] <- roots((obj$multi.model)[[j]])
        if (mod) r[[j]] <- Mod(r[[j]])
        print(r[[j]], digits=digits)
        cat("\n")
       }
    }
  invisible(r)
 }

############################################################################

#    functions for evaluating forecasts    <<<<<<<<<<<<<

# Class "horizon.forecasts"  has 
#  multiple horizons forecasts calculated from every data point.  
#  This is primarily used for calculating forecast errors at different horizons.

############################################################################

#    methods for horizon.forecasts        <<<<<<<<<<<<<

############################################################################

is.horizon.forecasts <- function(obj) { inherits(obj,"horizon.forecasts") }

test.equal.horizon.forecasts <- function(obj1, obj2, fuzz=1e-14)
{# N.B. models are not compared (so equivalent models will compare true)
 r <- all(dseclass(obj1) == dseclass(obj2))
 if (r) r <- test.equal(obj1$data, obj2$data)
 if (r) r <- all(obj1$horizons == obj2$horizons)
 if (r) r <- obj1$discard.before == obj2$discard.before
 if (r) r <- fuzz > max(abs(obj1$horizon.forecasts - obj2$horizon.forecasts))
 if (r) r <- test.equal(obj1$data, obj2$data)
 r
}

horizon.forecasts <- function(model, data, ...) UseMethod("horizon.forecasts")

horizon.forecasts.TSestModel <- function (model, data=NULL, ...) 
{horizon.forecasts(TSmodel(model), if (is.null(data)) TSdata(model) else data, ...)
}


horizon.forecasts.TSdata <- function(data,model, ...)
{ horizon.forecasts.TSmodel(model, data, ...)}

horizon.forecasts.TSmodel <- function(model, data, horizons=1:4,
	 discard.before=minimum.startup.lag(model), compiled=.DSECOMPILED)
{# calculate multiple "horizon"-step ahead forecasts 
 # ie. calculate forecasts but return only those indicated by horizons.
 #     Thus, for example, the result of
 #          horizon.forecasts(model, data horizons=c(1,5))    
 #     would be the one-step ahead and five step ahead forecasts.
 # The result is a list of class horizon.forecasts with elements model (a 
 #   TSmodel), data, horizons, discard.before, and horizon.forecasts.
 # horizon.forecasts is an array with three dimension: 
 #   c(length(horizons), dim(model$data)).
 # Projections are not calculated before discard.before or after
 #   the end of output.data(data).
 # Each horizon is aligned so that horizon.forecasts[h,t,] contains the forecast
 #   for the data point output.data(data)[t,] (from horizon[h] periods prior).
 
   if(!check.consistent.dimensions(model,data)) stop("dimension error\n")
   if (compiled) proj <- horizon.forecasts.compiled(model, data, 
                           horizons=horizons, discard.before=discard.before)
   else
     {TT <-periods(data)
      proj <- array(NA,c(length(horizons),dim(output.data(data))))
      for (t in discard.before:(TT-min(horizons)) )
        {horizons <- horizons[t <= (TT-horizons)]
         z <- l(model, data, sampleT=t, predictT=TT)$estimates$pred
         for (h in 1: length(horizons) )
             proj[h,t+horizons[h],] <- z[t+horizons[h],]
        }
     }
   dimnames(proj) <- list(NULL, NULL, output.series.names(data))
   proj <- list(model=model, data=data, horizons=horizons, 
                discard.before=discard.before, horizon.forecasts=proj)
   invisible(classed(proj, "horizon.forecasts" ))
}

# pbb2.6<- horizon.forecasts(l(bb2.vs.ls2[[3]]$model,aug10.KFmonitor.data), horizons=6, discard.before=)

# zf<- horizon.forecasts(l(bb2.vs.ls2[[3]]$model,aug10.KFmonitor.data), horizons=c(3,4), discard.before=220)

horizon.forecasts.compiled <- function(obj, ...) 
   UseMethod("horizon.forecasts.compiled")

horizon.forecasts.compiled.ARMA <- function( model, data, horizons=1:4,
	  discard.before=minimum.startup.lag(model))
{ if (discard.before < dim(model$A)[1] )
       warning(paste("Results may be spurious. discard.before should be set higher than the AR order (=", 
                   dim(model$A)[1]-1, ")."))
 horizons <- sort(horizons)
  p <- output.dimension(data)
  TT <- periods(data)
  proj <- array(0,c(length(horizons),TT,p))
  storage.mode(proj) <- "double"
  m <- input.dimension(model)
  if (m==0)
     {C <- array(0,c(1,p,1))    # can't pass 0 length array to compiled
      u <- matrix(0,TT,1)
     }
  else
     {C <-    model$C
      u <- input.data(data)
      if (discard.before < dim(C)[1] )
        warning(paste("Results may be spurious. discard.before should be set higher than the order of C (=", 
                      dim(C)[1]-1, ")."))
     }
  TREND <- model$TREND
  if (is.null(model$TREND)) TREND <- rep(0,p)
  is  <- max(m,p)
   .Fortran("rmaprj",
                  proj=proj,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  ey= as.double(array(0,dim(output.data(data)))), 
                  as.integer( m), 
                  as.integer( p) ,      
                  as.integer( dim(model$A)[1]),  # 1+order of A  
                  as.integer( dim(model$B)[1]),  # 1+order of B  
                  as.integer( dim(C)[1]),  # 1+order of C  
                  as.integer(TT),
                  as.double(u),
                  as.double(output.data(data)),         
                  as.double(model$A),  
                  as.double(model$B),   
                  as.double(C),
                  as.double(TREND),
                  as.integer(is),  # scratch array dim
                  as.double(matrix(0,is,is)),  # scratch array
                  as.double(matrix(0,is,is)),  # scratch array
                  as.double(rep(0,is))         # scratch array
             )$proj
}

horizon.forecasts.compiled.SS <- function( model, data, horizons=1:4,
	 discard.before=minimum.startup.lag(model))
{ horizons <- sort(horizons)
  p <- output.dimension(data)
  TT <- periods(data)
  proj <- array(0,c(length(horizons),TT,p))
  storage.mode(proj) <- "double"
     gain <- is.innov.SS(model)
     n <- dim(model$F)[2]
     if (discard.before <= n )
         warning(paste("discard.before should probably be set higher than the state dimension (=", n, ")."))
     if (is.null(model$G))
       {m<-0
        G<-matrix(0,n,1)       # can't call compiled with 0 length arrays
        u <- matrix(0,TT,1)
       }
     else
       {m <- dim(model$G)[2]
        G <-model$G
        u <- input.data(data)[1:periods(data),,drop=F]
       } 
     if (gain)     # K or Q,R can be NUll in model, which messes up compiled
       {K <-    model$K
        Q <-    matrix(0,1,1)      #not used
        R <-    matrix(0,1,1)      #not used
       }
     else
       {Q <-    model$Q
        if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square in compiled
        R <-    model$R
        K <-    matrix(0,n,p)      # this is used
       }
     if(is.null(model$z0)) z <-rep(0,n)   # initial state
     else  z <-model$z0
     if(is.null(model$P0)) P <- diag(1,n) # initial state tracking error 
     else  P <- model$P0              # this is not used in innov. models

     .Fortran("kfprj",
                  proj= proj, 
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  ey= as.double(matrix(0,TT,p)), 
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  as.integer(TT),  
                  as.double(u), 
                  as.double(output.data(data)),  
                  as.double(model$F),   
                  as.double(G),   
                  as.double(model$H),  
                  as.double(K), 
                  as.double(Q),      
                  as.double(R),    
                  as.logical(gain),
                  as.double(z),
                  as.double(P))$proj
}


forecasts.horizon.forecasts <- function(obj){obj$horizon.forecasts}

tfplot.horizon.forecasts <- function(x, start.=NULL, end.=NULL, select.series=NULL, names=output.series.names(x$data))
{#If select.series is not NULL then only indicated variables are plotted
 # if start. is null it is set to the beginning of the data.
 # if end. is null it is set to the end of the data.
   output <-output.data(x$data)
   if(is.null(names)) names <- rep(" ", dim(output)[2])
   if (is.null(select.series)) select.series <- 1: dim(output)[2]
   if (is.null(start.)) start. <- start(output)
   if (is.null(end.)) end. <- end(output)
   old.par <-par(mfcol = c(length(select.series), 1), mar= c(5.1,6.1,4.1,2.1))
   on.exit(par(old.par))
   tf <- tframe(output)
   for(i in select.series) 
     {#unclass below in because x$horizon.forecasts is not tframed and tbind
      #   complains if the frequencies do not match
      zz<- tframed(tbind(unclass(output)[,i],t((x$horizon.forecasts)[,,i])), tf)
      tfplot(tfwindow(zz,start=start.,end=end., warn=F), ylab =names[i])
      if(i == select.series[1]) title(main = "Actual data (solid)")
     }
   invisible()
}

 
############################################################################
#
#       methods for estimators.horizon.forecasts.wrt.data.   <<<<<<<<<<
#
############################################################################

estimators.horizon.forecasts <- function(data, 
                       estimation.sample=.5, horizons=1:12,quiet=F,
                       estimation.methods=NULL)
{ # estimation.sample indicates the part of the data to use for estimation.
  # If estimation.sample is less than or equal 1.0 it is
  # used to indicate the portion of points to use for estimation. 
  # Otherwise it should be an integer and is used to indicate the number
  # of points from the beginning of the sample to use for estimation. 
  
  if (is.null(estimation.methods)) stop("estimation.methods must be specified.")
  if (estimation.sample <= 1.0 )
     estimation.sample <- as.integer(round(estimation.sample*nrow(output.data(data))))
  r <- list(data=data, estimation.sample =estimation.sample, horizons=horizons,  
            estimation.methods=estimation.methods )

  r$multi.model <- estimate.models(data, estimation.sample=estimation.sample, 
              trend=F,quiet=quiet, 
	      estimation.methods=estimation.methods)$multi.model

  r$horizon.forecasts <- vector("list", length(estimation.methods))
  for (j in 1:length(estimation.methods))
    r$horizon.forecasts[[j]] <- horizon.forecasts(l(r$multi.model[[j]],data),
               horizons=horizons,
	       discard.before=minimum.startup.lag(r$multi.model[[j]]))
  classed(r, c("estimators.horizon.forecasts.wrt.data")) #? "horizon.forecasts")
}


############################################################################
#
#       methods for forecast.cov.   (including multiple models)<<<<<<<<<<
#
############################################################################

horizon.forecasts.forecast.cov <- function(obj,horizons=NULL, discard.before=NULL)
{# Calculate forecasts of an object for which cov has been calculated.
 # In a sense this is a step backward, but is sometimes useful to look at
 # forecasts after methods have been analysed on the basis of cov. 
 if(is.null(horizons))       horizons <- obj$horizons
 if(is.null(discard.before)) discard.before <- obj$discard.before
 if (!is.null(obj$model))
   {proj <- horizon.forecasts.TSmodel(obj$model, obj$data, horizons=horizons, 
                       discard.before=discard.before)
    dseclass(proj) <- "horizon.forecasts"
   }
 else if (!is.null(obj$multi.model))
   {proj <-vector("list", length(obj$multi.model))
    for (i in seq(length(obj$multi.model)))
      proj[[i]] <-horizon.forecasts.TSmodel(
             (obj$multi.model)[[i]], obj$data, 
             horizons=horizons, discard.before=discard.before)
    dseclass(proj) <- c("multi.model.horizon.forecasts","horizon.forecasts")
   }
 else  stop("Object does not include a model.\n")
 invisible(proj)
}

tfplot.multi.model.horizon.forecasts <- function(x, start.=NULL, end.=NULL, select.series=NULL)
 {for (i in seq(length(x)))
    {tfplot(x[[i]], start.=start., end.=end., select.series=select.series)
     cat("press return to continue>");key<-dsescan(what="");cat("\n")
    }
  invisible()
 }
# zz<-forecast.cov(zl, discard.before=1, horizons=1:12)
# z<-forecast.cov(l(mod3,simulate(mod3)), discard.before=20, horizons=1:12)
# zz<-forecast.cov(l(mod3,simulate(mod3)), discard.before=80, horizons=1:4)
# zzz<-forecast.cov(zz$model,zz$data, discard.before=zz$discard.before, horizons=zz$horizons)


forecast.cov <- function(obj, ...)
  {# Use model and data to construct the cov of predictions at horizons.
   # Discard predictions before (but not including) discard.before to remove 
   #    initial condition problems or do out-of-sample analysis.
   #  obj can be a TSestModel or 
   #        a TSmodel, in which case the second arg must be TSdata, or
   #          TSdata,  in which case the second arg must be a TSmodel.
   UseMethod("forecast.cov")
  }

forecast.cov.TSestModel <- function( obj, data=NULL, ...)
 {forecast.cov(TSmodel(obj), data=if(is.null(data)) TSdata(obj) else data, ...)}

forecast.cov.TSdata <- function( pred, data=NULL, horizons=1:12, discard.before=1, compiled=.DSECOMPILED)
{# Use pred$output as the predictions of data and calculate forecast.cov
 # This is mainly useful for a fixed prediction like zero or trend.
 # The calculation is dominated by sample effects: more points are
 #  dropped from the end for longer prediction horizons; the trend
 #  predictions are better for the first few periods.
 # With very large samples the result should be nearly constant for 
 # all horizons.
 # The default discard.before=1 should work ok for data, but is not 
 #    consistent with the value for model forecasts. When this routine is
 #    called by other functions the value will usually be overridden.
   horizons <- sort(horizons)
   p <- output.dimension(data)
   TT  <- periods(data)
   cov <- array(0,c(length(horizons), p,p))
   N <- rep(0,length(horizons))   # the sample size used at each horizon
   err <- pred$output - output.data(data)
   if (compiled)
     {storage.mode(cov) <-"double"
      storage.mode(err) <-"double"
      r <- .Fortran("datepr",
                  forecast.cov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  as.integer(p), 
                  predictT=as.integer(TT), 
                  as.double(err)) [c("forecast.cov","sample.size")]
     }
   else
     {for (t in discard.before:(TT-horizons[1]+1))
        {h <- t-1+horizons[(t-1+horizons) <= TT]
         e <- err[h,,drop=F]
         for (k in 1:length(h))
           {N[k] <- N[k]+1
            cov[k,,] <- cov[k,,]*((N[k]-1)/N[k]) + e[k,] %o% e[k,]/N[k] 
           }
        }
       r <- list( forecast.cov=cov, sample.size=N)
     }
  dimnames(r$forecast.cov) <- list(paste("horizon",as.character(horizons)),NULL,NULL)
  r$forecast.cov <- list(r$forecast.cov)
  r <- append(r, list(pred=pred, data=data, model=NULL, horizons=horizons, 
                      discard.before=discard.before))
  classed(r, "forecast.cov")
}

TSmodel.forecast.cov <- function(obj, select=1)
  {if (is.null(obj$multi.model)) NULL else obj$multi.model[[select]]}

TSdata.forecast.cov <- function(obj) {obj$data}

forecast.cov.TSmodel <- function(obj, ..., data=NULL, discard.before=NULL,
       horizons=1:12, zero=F, trend=F, estimation.sample= periods(data),
       compiled=.DSECOMPILED)
{# Calculate the forecast cov of models in list(obj, ...) with data.
 # Using obj, ... instead of just something like model.list make argument
 # matching a bit messier, but means the method gets called for a single 
 #  TSmodel (obj), or for a list of TSmodels (obj, ...), without making the 
 #  list into a class, etc.
 # This is just multiple applications of  forecast.cov.single.TSmodel
 # discard.before is an integer indicating the number of points in the
 #   beginning of forecasts to discard before calculating covariances.
 #   If it is the default, NULL, then the default (minimum.startup.lag) will
 #   be used for each model and the default (1) will be used for trend and zero.
 # If zero  is T then forecast.cov is also calculated for a forecast of zero.
 # If trend is T then forecast.cov is also calculated for a forecast of a linear
 #   trend using data to estimation.sample.
  if (is.null(data)) stop("data= must be supplied.")
  model.list <- list(obj, ...)
  r <- list(data=data, horizons=horizons, discard.before =discard.before)
  if (is.TSmodel(model.list)) model.list <- list(model.list)
  r$forecast.cov <-vector("list", length(model.list))
  i <-0  
  for (model in model.list)
      {i <- i+1
       if (is.null(discard.before))
             rn <-  forecast.cov.single.TSmodel(TSmodel(model), data, 
                           horizons=horizons, compiled=compiled)
       else  rn <-  forecast.cov.single.TSmodel(TSmodel(model), data, 
                           horizons=horizons, discard.before=discard.before, 
                           compiled=compiled)
       #  $ in the following causes problems for some reason
       r$forecast.cov[[i]] <- rn$forecast.cov
       r$sample.size   <- rn$sample.size
      }
  if (trend)
     {y <- output.data(data)[1:estimation.sample,]
      pred <- cbind(1,1:periods(data)) %*%
                              (lsfit(1:estimation.sample, y)$coef)
      if (is.null(discard.before))
         r$forecast.cov.trend <- (forecast.cov.TSdata(list(output=pred), data,
             horizons=horizons)$forecast.cov)[[1]]
      else
         r$forecast.cov.trend <- (forecast.cov.TSdata(list(output=pred), data,
           horizons=horizons,discard.before=discard.before)$forecast.cov)[[1]]
     }
  if (zero)
     {if (is.null(discard.before))
        r$forecast.cov.zero <- (forecast.cov.TSdata(
             list(output=array(0,dim(output.data(data)))), data, 
             horizons=horizons)$forecast.cov)[[1]]
      else
        r$forecast.cov.zero <- (forecast.cov.TSdata(
             list(output=array(0,dim(output.data(data)))), data,
           horizons=horizons,discard.before=discard.before)$forecast.cov)[[1]]
     }
  r$multi.model <- model.list
  classed(r, c("forecast.cov.wrt.data", "forecast.cov"))
}

is.forecast.cov.wrt.data <- function(obj){inherits(obj,"forecast.cov.wrt.data")}

forecast.cov.single.TSmodel <- function( model, data=NULL, horizons=1:12, 
          discard.before=minimum.startup.lag(model), compiled=.DSECOMPILED)
{ if(!check.consistent.dimensions(model,data)) stop("dimension error.")
  if (discard.before < 1) stop("discard.before cannot be less than 1.")
  horizons <- sort(horizons)
  names <- series.names(data)$output
  if (compiled) 
     r <- forecast.cov.compiled(model, data, horizons=horizons,
     		 discard.before=discard.before)
  else
    { p <- output.dimension(data)
      shf <- start.shift(model,data) #,y0=y0)
      TT  <-periods(data)-(shf$shift)*(shf$lags+shf$terminal.periods)
      cov <- array(0,c(length(horizons), p,p))
      N <- rep(0,length(horizons))   # the sample size used at each horizon
  # there is a problem here with troll models trying to simulate further than
  #  the database allows. (after many steps.)
      for (t in discard.before:(TT-horizons[1]+1))
        {pred <- l(model, data, sampleT=t, predictT=TT, result="pred")
         # Eliminate longer horizons as data runs out.
         # This assumes HORIZ is sorted in ascending order.
         h <- t-1+horizons[(t-1+horizons) <= TT]
         e <- pred[h,,drop=F]- output.data(data)[h,,drop=F]
         for (k in 1:length(h))
           {N[k] <- N[k]+1
            cov[k,,] <- cov[k,,]*((N[k]-1)/N[k]) + e[k,] %o% e[k,]/N[k] 
           }
        }
       r <- list( forecast.cov=cov, sample.size=N)
     }
  dimnames(r$forecast.cov) <- list(paste("horizon",as.character(horizons)),names,names)
#  old:
# The following puts the cov in a sub list. This seems unnecessary for a single
#   cov, but means the same structure can be used with multiple model covs.
#  r$forecast.cov <- list(r$forecast.cov)
#  r <- append(r, list(model=model, data=data, horizons=horizons, 
#                     discard.before=discard.before))
#  class(r) <- "forecast.cov"
 r
}

forecast.cov.compiled <- function(obj, ...) 
   {if (!exists(paste("forecast.cov.compiled.", dseclass(obj)[1], sep="")))
stop("compiled code for this model class is not available. Try forecast.cov( ..., compiled=F)")
    UseMethod("forecast.cov.compiled")
   }

forecast.cov.compiled.ARMA <- function( model, data, horizons=1:12 , discard.before=minimum.startup.lag(model))
{ if (discard.before < dim(model$A)[1] )
       warning(paste("Results may be spurious. discard.before should be set higher than the AR order (=",
                    dim(model$A)[1]-1, ")."))
  horizons <- sort(horizons)
  p <- output.dimension(data)
  TT <- periods(data)
  cov <- array(0,c(length(horizons), p,p))
  N <- rep(0,length(horizons))   # the sample size used of each horizon
  m <- dim(model$C)[3]
  if (is.null(model$C))
     {m <- 0
      C <- array(0,c(1,p,1))    # can't pass 0 length array to compiled
      u <- matrix(0,TT,1)
     }
  else
     {C <-    model$C
      m <- dim(model$C)[3]
      u <- input.data(data)
      if (discard.before < dim(C)[1] )
        warning(paste("Results may be spurious. discard.before should be set higher than the order of C (=", dim(C)[1]-1, ")."))
     }
  TREND <- model$TREND
  if (is.null(model$TREND)) TREND <- rep(0,p)
  storage.mode(cov) <-"double"
  is  <- max(m,p)
  .Fortran("rmaepr",
                  forecast.cov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  pred= as.double(array(0,dim(output.data(data)))), 
                  as.integer( m), 
                  as.integer( p) ,      
                  as.integer( dim(model$A)[1]),  # 1+order of A  
                  as.integer( dim(model$B)[1]),  # 1+order of B  
                  as.integer( dim(C)[1]),  # 1+order of C  
                  predictT=as.integer(TT),
                  as.integer(nrow(output.data(data))), 
                  as.double(u),
                  as.double(output.data(data)),         
                  as.double(model$A),  
                  as.double(model$B),   
                  as.double(C),
                  as.double(TREND),
                  as.integer(is),  # scratch array dim
                  as.double(matrix(0,is,is)),  # scratch array
                  as.double(matrix(0,is,is)),  # scratch array
                  as.double(rep(0,is))         # scratch array
              )[c("forecast.cov","sample.size")]
}

forecast.cov.compiled.innov <- function(obj, ...)
  {forecast.cov.compiled.SS(obj, ...)}
forecast.cov.compiled.non.innov <- function(obj, ...) 
  {forecast.cov.compiled.SS(obj, ...)}

forecast.cov.compiled.SS <- function( model, data, horizons=1:12 , discard.before=minimum.startup.lag(model))
{ horizons <- sort(horizons)
  p <- output.dimension(data)
  TT <- periods(data)
  cov <- array(0,c(length(horizons), p,p))
  N <- rep(0,length(horizons))   # the sample size used at each horizon
     gain <- is.innov.SS(model)
     n <- dim(model$F)[2]
     if (discard.before <= n )
       warning(paste("discard.before should probably be set higher than the state dimension (=", n, ")."))
     if (is.null(model$G))
       {m<-0
        G<-matrix(0,n,1)       # can't call compiled with 0 length arrays
        u <- matrix(0,TT,1)
       }
     else
       {m <- dim(model$G)[2]
        G <-model$G
        u <- input.data(data)
       } 
     if (gain)     # K or Q,R can be NUll in model, which messes up compiled
       {K <-    model$K
        Q <-    matrix(0,1,1)      #not used
        R <-    matrix(0,1,1)      #not used
       }
     else
       {Q <-    model$Q
        if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square in compiled
        R <-    model$R
        K <-    matrix(0,n,p)      # this is used
       }
     if(is.null(model$z0)) z <-rep(0,n)   # initial state
     else  z <-model$z0
     if(is.null(model$P0)) P <- diag(1,n) # initial state tracking error 
     else  P <- model$P0              # this is not used in innov. models

     storage.mode(cov) <-"double"
     r <- .Fortran("kfepr",
                  forecast.cov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  pred= as.double(array(0,dim(output.data(data)))), 
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  predictT=as.integer(TT), 
                  as.integer(nrow(output.data(data))),  
                  as.double(u), 
                  as.double(output.data(data)),  
                  as.double(model$F),   
                  as.double(G),   
                  as.double(model$H),  
                  as.double(K), 
                  as.double(Q),      
                  as.double(R),    
                  as.logical(gain),
                  as.double(z),
                  as.double(P)) [c("forecast.cov","sample.size")]
  r
}

is.forecast.cov <- function(obj)
{ inherits(obj,"forecast.cov") }

print.forecast.cov <- function(x, digits=4)
{for (i in 1:dim((x$forecast.cov)[[1]])[3]) 
     {cat("   ",dimnames(x$forecast.cov)[[1]][i], "\n")
      z <- NULL
      for (j in 1:length(x$forecast.cov) )
         z <- tbind(z,  (x$forecast.cov)[[j]][,i,i])
      print(z, digits=digits)
     }
 invisible(x)
}


summary.forecast.cov <- function(object, horizons=object$horizons, 
    select.series=seq(output.dimension(object$data)))
{ names <- output.series.names(object$data)
 if(!is.numeric(select.series)) select.series <- match(select.series, names)
 names <- names[select.series]
 descriptions  <-  vector(NA,    length(object$multi.model))
 summary.stats <-  vector("list",length(object$multi.model))

 for (i in seq(length(summary.stats)))
   {descriptions[i] <- object$multi.model[[i]]$description
    z <- NULL
    for (h in seq(length(horizons))) z <- rbind(z,
             diag(object$forecast.cov[[i]][h,select.series,select.series])^0.5)
    dimnames(z) <- list(paste("S.D.horizon", horizons), names)
    summary.stats[[i]] <- z
   }
  classed(list(  # constructor summary.forecast.cov
     class=dseclass(object),
     horizons=length(object$horizons),
     models=length(object$multi.model),
     output.series.names=output.series.names(object$data),
     descriptions=descriptions,
     names=names,
     summary.stats=summary.stats,
     nxt=NextMethod("summary")),
  "summary.forecast.cov")
 }



print.summary.forecast.cov <- function(x, digits=options()$digits)
 {cat("class: ", x$class,"   ")
  cat( length(x$horizons), " horizons\n")
  cat(length(x$models),"models\n")
  for (i in seq(length(x$models)))
    {cat("Model", i, x$description[i],  "\n")
     cat("   variable", x$names,"\n")
     print(summary.stats[[i]], digits = digits)   
    }
  cat("\n")
  print(x$nxt)
  invisible(x)
 }





test.equal.forecast.cov <- function(obj1, obj2, fuzz=1e-14)
{if (is.null(obj1$rng)) ok <- T
 else ok <- test.equal(obj1$rng , obj2$rng)
 if (ok & !is.null(obj1$forecast.cov.true) )
  {if (is.null(obj2$forecast.cov.true)) ok <-F
   ok <- fuzz > max(abs(obj1$forecast.cov.true-obj2$forecast.cov.true))
  }
 if (ok & !is.null(obj1$forecast.cov.zero)) 
  {if (is.null(obj2$forecast.cov.zero)) ok <-F
   else ok <- fuzz > max(abs(obj1$forecast.cov.zero-obj2$forecast.cov.zero))
  }
 if (ok & !is.null(obj1$forecast.cov.trend)) 
  {if (is.null(obj2$forecast.cov.trend)) ok <-F
   else ok <- fuzz > max(abs(obj1$forecast.cov.trend-obj2$forecast.cov.trend))
  }
 for (i in 1:length(obj1$forecast.cov))
   {if (ok & !is.null((obj1$forecast.cov)[[i]])) 
         {if (is.null((obj2$forecast.cov)[[i]])) ok <-F 
          else ok <- fuzz > 
               max(abs((obj1$forecast.cov)[[i]]-(obj2$forecast.cov)[[i]]))
      }
   }
 ok
}


total.forecast.cov <- function(obj, select=NULL)
{if (is.null(select)) select <-1:dim((obj$forecast.cov)[[1]])[2]
 N <- c( dim((obj$forecast.cov)[[1]])[1] ,1,1)
 for (j in 1:length(obj$forecast.cov) )
   {z <- apply((obj$forecast.cov)[[j]],1,diag)
    # $ causes problems
    obj$forecast.cov[[j]] <- array(apply(z[select,],2,sum), N)
   }
 if(!is.null(obj$forecast.cov.true))
   {z <- apply(obj$forecast.cov.true,1,diag)
    obj$forecast.cov.true <- array(apply(z[select,],2,sum), N)
   }
 if(!is.null(obj$forecast.cov.zero))
   {z <- apply(obj$forecast.cov.zero,1,diag)
    obj$forecast.cov.zero <- array(apply(z[select,],2,sum), N)
   }
 if(!is.null(obj$forecast.cov.trend)) 
   {z <- apply(obj$forecast.cov.trend,1,diag)
    obj$forecast.cov.trend <- array(apply(z[select,],2,sum), N)
   }
 invisible(obj)
}



tfplot.forecast.cov <- function(x, select.series = 1:dim(x$forecast.cov[[1]])[2], 
    select.cov = 1:length(x$forecast.cov), select.true = T, 
    select.zero = T, select.trend = T, y.limit = NULL, line.labels = F, 
    lty = NULL, Legend = NULL, Title = NULL, graphs.per.page = 5, 
    ...) 
{
    # ... should be arguments to par().
    # select.cov indicates which covariances to display 
    #  (ie. which model or estimation method)
    # select.series, if specified, indicates which series to display.
    #  cex= can be passed as an argument to change character print size.
    #  If lty is NULL (default) it is set to
    #    seq(length(select.cov) +select.true+select.zero+select.trend),
    #  and corrected if these are T but not in the object.
    p <- dim((x$forecast.cov)[[1]])[2]
    Ngraph <- 1 + min(length(select.series), graphs.per.page)
    old.par <- par(mfcol = c(Ngraph, 1), mar = c(5.1, 6.1, 4.1, 
        2.1))
    on.exit(par(old.par))
    par(...)
    if (is.null(lty)) 
        lty <- seq(length(select.cov) + 
            (select.true  & !is.null(x$forecast.cov.true)) + 
            (select.zero  & !is.null(x$forecast.cov.zero)) + 
            (select.trend & !is.null(x$forecast.cov.trend)))
    names <- dimnames((x$forecast.cov)[[1]])[[2]]
    if (is.null(names)) 
        names <- paste("variable", 1:p)
    for (i in select.series) {
        z <- matrix(0, length(x$horizons), length(select.cov))
        for (j in 1:length(select.cov))
            z[, j] <- (x$forecast.cov)[[select.cov[j]]][, i, i]
        if (select.trend & !is.null(x$forecast.cov.trend)) 
            z <- tbind((x$forecast.cov.trend)[, i, i], z)
        if (select.zero & !is.null(x$forecast.cov.zero)) 
            z <- tbind((x$forecast.cov.zero)[, i, i], z)
        if (select.true & !is.null(x$forecast.cov.true)) 
            z <- tbind((x$forecast.cov.true)[, i, i], z)
        show <- 1:length(select.cov)
        if (!is.null(y.limit)) {
            z[z > y.limit] <- NA
            #Legend may be messed up for multi-v. case
            show <- !apply(is.na(z), 2, all)
            z <- z[, show]
        }
        if (is.character(i)) 
            ylab <- i
        else ylab <- names[i]
        #cex, mar set by par above
        matplot(x$horizons, z, type = "l", lty = lty, xlab = "horizon", 
            ylab = ylab)
        if (line.labels) {
            labels <- select.cov[show]
            if (!is.null(x$selection.index)) 
                labels <- x$selection.index[labels]
            text(dim(z)[1], z[dim(z)[1], ], labels)
        }
        if ((i == 1) & (!is.null(Title))) 
            title(main = Title)
    }
    if (is.null(Legend)) {
        Legend <- paste("prediction covariance", select.cov[show])
        if (!is.null(x$variable.index)) 
           {nm <- c(output.series.names(x$all.data), 
                  input.series.names(x$all.data))
            for (i in 1:length(select.cov))
                Legend[i] <- paste(Legend[i], "using", 
                   paste(nm[x$variable.index[select.cov[i],]], collapse=" "))
           }
        if (!is.null(x$selection.index)) 
            Legend <- paste(Legend, "=", x$selection.index[show])
        if (select.trend & !is.null(x$forecast.cov.trend)) 
            Legend <- c("trend", Legend)
        if (select.zero & !is.null(x$forecast.cov.zero)) 
            Legend <- c("zero", Legend)
        if (select.true & !is.null(x$forecast.cov.true)) 
            Legend <- c("true", Legend)
    }
    par(mfg = c(Ngraph, 1, Ngraph, 1)) 
    if (is.R()) 
       {box(col = 0) #  temp Rbug workaround
        legend((par()$usr)[1:2], (par()$usr)[3:4], Legend, lty = lty, 
            col = 1:6, bty = "y")
       }
    else legend((par()$usr)[1:2], (par()$usr)[3:4], Legend, lty = lty, 
               col = lty, bty = "y") # this is a bit of an S/R comp. issue.
        #bty is box
    invisible()
}


############################################################################
#
#       end
#
############################################################################

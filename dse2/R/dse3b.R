
############################################################################
#
#   methods for "forecast.cov.estimators.wrt.data", "forecast.cov"  <<<<<<<<<<
#                  (multiple estimators,  given data)
#
############################################################################
#z <-out.of.sample.forecast.cov.estimators.wrt.data(zl$data,
#       estimation.methods = list(est.VARX.ar=list(max.lag=2),
#                                 est.VARX.ls=list(max.lag=2)))

#z <-out.of.sample.forecast.cov.estimators.wrt.data(data,
#             estimation.sample=.5,trend=T, zero=T,
#             estimation.methods = list(
#                   est.VARX.ls=list(max.lag=4),
#                   est.wt.variables=list(
#                       variable.weights=c(1,1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5),
#                       estimation.methods=list(est.VARX.ls=list(max.lag=4))))


out.of.sample.forecast.cov.estimators.wrt.data <- function(data, zero=F, trend=F,
                       estimation.sample=.5, horizons=1:12,quiet=F,
                       estimation.methods=NULL, compiled=.DSECOMPILED)
{ # estimation.sample indicates the portion of the data to use for estimation.
  #If estimation.sample is an integer then it is used to indicate the number
  # of points in the sample to use for estimation. If it is a fracton it is
  # used to indicate the portion of points to use for estimation. The remainder
  # of the sample is used for evaluating forecasts.
  
  if (estimation.sample < 1.0 )
     estimation.sample <- as.integer(round(estimation.sample*nrow(output.data(data))))
  discard.before <- 1+estimation.sample
  forecast.cov.estimators.wrt.data(data, estimation.sample, discard.before,
                       horizons=horizons, zero=zero, trend=trend, quiet=quiet,
                       estimation.methods=estimation.methods, compiled=compiled)
}


forecast.cov.estimators.wrt.data <- function(data, estimation.sample=NULL, 
                       compiled=.DSECOMPILED, discard.before=10,
                       horizons=1:12, zero=F, trend=F,quiet=F,
                       estimation.methods=NULL)
{# Calculate the forecasts cov of models estimated from data with estimation
 #   methods indicated by estimation.methods  (see estimate.models).
 # estimation.sample is an integer indicating the number of points in the
 #     sample to use for estimation. If it is NULL the whole sample is used.
 # discard.before is an integer indicating 1+the number of points in the
 #     beginning of forecasts to discard for calculating covariances.
 # If zero  is T then forecast.cov is also calculated for a forecast of zero.
 # If trend is T then forecast.cov is also calculated for a forecast of a linear trend.

  r <- list(data=data, estimation.sample =estimation.sample,
            horizons=horizons, discard.before =discard.before, 
            estimation.methods=estimation.methods)
  models <-estimate.models(data, estimation.sample=estimation.sample, 
                       trend=trend,quiet=quiet,
                       estimation.methods=estimation.methods)
  r$multi.model <- models$multi.model
  if (!is.null(estimation.methods))
    {r$forecast.cov <- vector("list", length(estimation.methods))
     for (j in 1:length(estimation.methods))
        {rn <-  forecast.cov.single.TSmodel(r$multi.model[[j]], data, 
             compiled=compiled, discard.before=discard.before, horizons=horizons)
         r$forecast.cov[[j]] <- rn$forecast.cov
         r$sample.size   <- rn$sample.size
        }
    }
  if (zero)
     {r$forecast.cov.zero <-forecast.cov.TSdata(
             list(output=array(0,dim(output.data(data)))), data, 
                discard.before=discard.before, 
                horizons=horizons)$forecast.cov[[1]]
     }
  if (trend)
     {pred <- cbind(1,1:periods(data)) %*% models$trend.coef
      r$forecast.cov.trend <- forecast.cov.TSdata(list(output=pred), data, 
              discard.before=discard.before, 
              horizons=horizons)$forecast.cov[[1]]
     }
  classed(r, c("forecast.cov.estimators.wrt.data", "forecast.cov")) # constructor
}

is.forecast.cov.estimators.wrt.data <- function(obj)
  {inherits(obj,"forecast.cov.estimators.wrt.data")}

combine.forecast.cov.estimators.wrt.data <- function(e1,e2)
  {if(! test.equal(e1$data, e2$data)) 
       warning("data is not the same. Second set suppressed.")
   if(! all(e1$estimation.sample == e2$estimation.sample)) 
       warning("estimation.sample's are not the same. Second one suppressed.")
   if(! all(e1$horizon == e2$horizon)) 
       stop("horizon's are not the same.")
   if(e1$discard.before != e2$discard.before) 
       warning("discard.before's are not the same. Second one suppressed.")
   e1$forecast.cov <- append(e1$forecast.cov, e2$forecast.cov)
   e1$estimation.methods <- append(e1$estimation.methods, e2$estimation.methods)
# fix   e1$multi.model <- append(e1$multi.model, e2$multi.model)  
   e1
}


extract.forecast.cov.estimators.wrt.data <- function(e,n)
  {UseMethod("extract.forecast.cov")}
  
extract.forecast.cov.estimators.wrt.data <- function(e,n)
  {# select indicated forecast.cov
   e$forecast.cov <- e$forecast.cov[[n]]
   e$estimation.methods <- e$estimation.methods[[n]]
   e$multi.model        <- e$multi.model[[n]]  
   e
}

tfplot.forecast.cov.estimators.wrt.data <- function(x, 
    select.series=1:dim(x$forecast.cov[[1]])[2], 
    select.cov=1:length(x$forecast.cov),
    select.zero=T, select.trend=T,
    lty=NULL,  ...)  # ,  lty=1:5
{# ... should be arguments to par(). See tfplot.forecast.cov for more details.
Legend<- paste(names(x$estimation.methods), x$estimation.methods)[select.cov]
 if(select.trend & !is.null(x$forecast.cov.trend))
       Legend  <- c("trend",Legend)
 if(select.zero  & !is.null(x$forecast.cov.zero))
       Legend  <- c( "zero",Legend)
 tfplot.forecast.cov(x, select.series=select.series, lty=lty,
        select.cov=select.cov, select.true=FALSE,
        select.zero=select.zero, select.trend=select.trend, Legend=Legend, 
        Title="Prediction variance relative to given data.",
        ...)
 invisible()
}

#graph.forecast.cov.estimators.wrt.data <- function( ..., select=NULL)
#{obj <- list( data = list(...)[[1]]$data)
# i <-0
# for (obji in list(...) )
#   {i <- i+1
#    obji$data <-NULL
#    na <- paste("obj$obj",as.character(i), sep="")
#    na <- paste(paste(na, names(obji), sep=""), " <-obji[[j]]" )
#    for (j in 1:length(obji)) eval(parse(text=na[j]))
#   }
# class(obj) <- "forecast.cov.estimators.wrt.data"
# tfplot.forecast.cov.estimators.wrt.data(obj, select=select)
# invisible()
#}

#date.ts <- function(x,i)
#{# date of ith position in time series x
# s <- start(x)
# tsp.x <-tsp(x)
# p <-s[2]+i
# y <- s[1]+ ((p-1) %/% tsp.x[3])
# p <- 1+(p %% tsp.x[3])
# c(y,p)
#}


############################################################################
#
#   methods for "forecast.cov.wrt.true", "forecast.cov"  <<<<<<<<<<
#    given true model, evaluate multiple estimation techniques
#    with multiple simulations for estimation and
#         multiple simulations for forecast
#
############################################################################

forecast.cov.wrt.true <- function( models, true.model, 
        pred.replications=1, simulation.args=NULL, quiet=F, 
        rng=NULL, Spawn=.SPAWN, compiled=.DSECOMPILED,
        horizons=1:12, discard.before=10, trend=NULL, zero=NULL)
{# models should be a list of models
 # The true model is used to generate more
 # data and for each generated data set the forecasts of the 
 # models are evaluated against the simulated data.
 # if trend is not null it is treated as a model output (forecast) and
 # should be the same dimension as a simulation of the models with 
 # simulation.args. If zero is not null a zero forecast is also evaluated.

 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
      
 if (Spawn & (pred.replications > 1))
   {if(!quiet) 
      cat("Spawning processes to calculate ", pred.replications,
            " forecast replications.\n")
    rep.forloop <- function(models, true.model, simulation.args,
                         horizons, discard.before, zero, trend, compiled=.DSECOMPILED)
      {data<-do.call("simulate",append(list(true.model), simulation.args))
       r <- NULL
       for (j in 1:length(models))
              {r <- c(r, forecast.cov.single.TSmodel(models[[j]],data,
                                  compiled=compiled, horizons=horizons, 
                                  discard.before=discard.before)$forecast.cov)
              }
         r.true <- forecast.cov.single.TSmodel(true.model,data, compiled=compiled,
                               horizons=horizons,
                               discard.before=discard.before)$forecast.cov
         if (is.null(trend)) r.trend <- NULL
         else  r.trend <-forecast.cov.TSdata(list(output=trend),
            data, discard.before=discard.before,horizons=horizons)$forecast.cov[[1]]
         if(is.null(zero)) r.zero <- NULL
         else r.zero <- forecast.cov.TSdata(
              list(output=array(0,dim(output.data(data)))), data, 
                discard.before=discard.before, horizons=horizons)$forecast.cov[[1]]
         c(dim(r.true),r.true, r.zero, r.trend,r)
       }

    assign("rep.forloop", rep.forloop, where = 1)
    assign("rep.forloop.n", pred.replications, where = 1)
    assign("rep.forloop.result", 0, where = 1)
    assign("rep.forloop.true.model", true.model, where = 1)
    assign("rep.forloop.simulation.args", simulation.args, where = 1)
    assign("rep.forloop.models", models, where = 1)
    assign("rep.forloop.horizons", horizons, where = 1)
    assign("rep.forloop.discard.before", discard.before, where = 1)
    assign("rep.forloop.trend", trend, where = 1)
    assign("rep.forloop.zero", zero, where = 1)
    assign("rep.forloop.compiled", compiled, where = 1)
    on.exit(remove(c("rep.forloop", "rep.forloop.i", "rep.forloop.n",
       "rep.forloop.models",  "rep.forloop.true.model",
       "rep.forloop.simulation.args","rep.forloop.result",
       "rep.forloop.horizons", "rep.forloop.discard.before", 
       "rep.forloop.trend","rep.forloop.zero","rep.forloop.compiled"),where = 1))

    For(rep.forloop.i = 1:rep.forloop.n, 
       rep.forloop.result <- rep.forloop.result +
          rep.forloop(rep.forloop.models, rep.forloop.true.model,
          rep.forloop.simulation.args, rep.forloop.horizons,
          rep.forloop.discard.before, rep.forloop.zero, rep.forloop.trend,
          rep.forloop.compiled),
       first=options(warn=-1), sync = T)

    names <- list(paste("horizon",as.character(horizons)),NULL,NULL)
    result  <- rep.forloop.result/pred.replications
    d <- result[1:3]  # this is not a very elegant way to pass this info.
    l <- prod(d)
    r.true  <- array(result[4:(3+l)], d)
    dimnames(r.true)  <- names
    lj <- 1
    if (is.null(zero)) r.zero <- NULL
    else
      {r.zero  <- array(result[(4+l):(3+2*l)], d)
       lj <-lj+1
       dimnames(r.zero)  <- names
      }
    if (is.null(trend)) r.trend <- NULL
    else
      {r.trend  <- array(result[(4+lj*l):(3+(1+lj)*l)], d)
       lj <-lj+1
       dimnames(r.trend) <- names
      }
    r <- vector("list", length(models))
    for (j in 1:length(models))
      {r[[j]] <- array(result[(4+lj*l):(3+(1+lj)*l)], d) 
       lj <-lj+1
       dimnames(r[[j]])  <- names
      }
     }
   else
     {r <-vector("list",length(models))
      r.true <- 0
      if (!is.null(zero)) r.zero <- 0
      else r.zero <- NULL
      if (!is.null(trend)) r.trend<- 0
      else r.trend <- NULL
      for (i in 1:pred.replications)
        {data<-do.call("simulate",append(list(true.model), simulation.args))
         for (j in 1:length(models))
              {fc <- forecast.cov.single.TSmodel(models[[j]],data,
                                  compiled=compiled, horizons=horizons, 
                                  discard.before=discard.before)$forecast.cov
               if (i == 1) r[[j]] <- fc
               else        r[[j]] <- r[[j]] + fc
              }
         r.true <- r.true+forecast.cov.single.TSmodel(TSmodel(true.model),data, 
                               compiled=compiled, horizons=horizons,
                               discard.before=discard.before)$forecast.cov
         if (!is.null(trend))
           r.trend <-r.trend+forecast.cov.TSdata(list(output=trend), data,
                                       discard.before=discard.before,
                                       horizons=horizons)$forecast.cov[[1]]
         if(!is.null(zero))
            r.zero <- r.zero + forecast.cov.TSdata(
                list(output=array(0,dim(output.data(data)))), data, 
                 discard.before=discard.before, horizons=horizons)$forecast.cov[[1]]
        }
      for (j in 1:length(models))  r[[j]] <- r[[j]]/pred.replications
      r.true  <-  r.true/ pred.replications
      if (!is.null(zero)) r.zero  <-  r.zero/ pred.replications
      if (!is.null(trend)) r.trend <-  r.trend/pred.replications
     }
   invisible(classed(  # constructor (forecast.cov.wrt.true)
         list(forecast.cov=r, forecast.cov.true=r.true, 
           forecast.cov.zero=r.zero, forecast.cov.trend=r.trend,
           multi.model=models,
           rng=rng, version=version,
           pred.replications=pred.replications,
           horizons=horizons, discard.before=discard.before),
        c("forecast.cov.wrt.true", "forecast.cov")))
}


forecast.cov.estimators.wrt.true <- function(true.model, Spawn=.SPAWN, rng=NULL,
                       simulation.args=NULL,
                       est.replications=2, pred.replications=2,
                       discard.before=10, horizons=1:12,quiet=F,
                       estimation.methods=NULL, compiled=.DSECOMPILED)
{# Calculate the forecasts cov of models estimated from simulations of 
 # true.model with estimation methods indicated by estimation.methods (see 
 #       estimate.models). 
 # discard.before is an integer indicating 1+the number of points in the
 # beginning of forecasts to discard for calculating forecast covariances.
 # The returned results has element
 #  $forecast.cov.true  $forecast.cov.zero $forecast.cov.trend containing 
 #    covariances averaged over estimation replications and simulation
 #    replications (forecasts will not change but simulated data will).
 #  $forecast.cov a list of the same length as estimation.methods with each
 #    element containing covariances averaged over estimation replications 
 #    and simulation replications.
 #  $estimated.models a list of length est.replications, with each elements as
 #    returned by estimate.models, thus each element has $multi.model as a
 #    subelement containing models for different estimation techniques.  
 #    So, eg.   $estimated.models[[2]]$multi.model[[1]]  in the result will
 #    be the model from the first estimation technique in the second replication. 

 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 
 estimated.models <- vector("list", est.replications)
 for (i in 1:est.replications)
        {data<-do.call("simulate",append(list(true.model), simulation.args))
         models <-estimate.models(data, trend=T,quiet=quiet,
                       estimation.methods=estimation.methods)
         estimated.models[[i]] <- models
         rn <- forecast.cov.wrt.true( models$multi.model, true.model, 
                    pred.replications=pred.replications, zero=T, quiet=quiet,
                    simulation.args=simulation.args, Spawn=Spawn,
                    horizons=horizons, discard.before=discard.before,
                    trend=cbind(1,1:periods(data)) %*% models$trend.coef,
                    compiled=compiled)
         if (i==1)
           r<-rn[c("forecast.cov","forecast.cov.true",
                   "forecast.cov.zero","forecast.cov.trend")]
         else
          {for (j in 1:length(estimation.methods))
              r$forecast.cov[[j]] <-  r$forecast.cov[[j]]*(i-1)/i + 
                                     rn$forecast.cov[[j]]/i
              r$forecast.cov.true  <- r$forecast.cov.true *(i-1)/i +  
                                     rn$forecast.cov.true/i
              r$forecast.cov.zero  <- r$forecast.cov.zero *(i-1)/i +  
                                     rn$forecast.cov.zero/i
              r$forecast.cov.trend <- r$forecast.cov.trend*(i-1)/i +  
                                     rn$forecast.cov.trend/i
          }
        }

  classed(append(r, # constructor forecast.cov.estimators.wrt.true
       list(true.model=true.model,estimation.methods=estimation.methods,
         estimated.models=estimated.models,
         rng=rng, version=version,
         horizons=horizons, 
         discard.before=discard.before, est.replications=est.replications,
         pred.replications=pred.replications, simulation.args=simulation.args)),
      c("forecast.cov.estimators.wrt.true", "forecast.cov") )
}

is.forecast.cov.estimators.wrt.true <- function(obj)
 {inherits(obj,"forecast.cov.estimators.wrt.true")}


print.forecast.cov.estimators.wrt.true <- function(x, digits=4)
{cat("forecast.cov.estimators.wrt.true\n")
 cat("essential data:", x$essential.data, "\n")
 cat("considering:", output.series.names(x$all.data), 
                      input.series.names(x$all.data), "\n")
 invisible(x)
}

summary.forecast.cov.estimators.wrt.true <- function(object, digits = 4)
 {conv <- list()
  Ms <- length(object$estimated.models)
  for (i in seq(Ms)) 
     conv<- append(conv,object$estimated.models[[i]]$multi.model[[1]]$converged)
  classed(list( # constructor summary.forecast.cov.estimators.wrt.true
    dseclass(object), 
    horizons=length(object$horizons), 
    Ms=Ms,
    conv=conv,
    nxt=NextMethod("summary")),
  "summary.forecast.cov.estimators.wrt.true")
 }


print.summary.forecast.cov.estimators.wrt.true <- function(x, 
      digits=options()$digits)
{cat("class: ", x[[1]], "   ")
 cat(x$horizons, " horizons\n")
 for (i in seq(x$Ms))
   {if (!is.null(x$conv))
      {cat("Estimated model", i)
       if(!x$conv) cat(" NOT")
       cat(" converged.\n")
   }  }
 print(x$nxt)
invisible(x)
}

roots.forecast.cov.estimators.wrt.true <- function(obj, digits=4, mod=F)
 {cat("Estimated models:\n")
  if (!is.null(obj$trend.coef)) cat("trend coef: ", obj$trend.coef, "\n")
  if (!is.null(obj$estimated.models))
    {r <- vector("list",length(obj$estimated.models)) 
     for (i in 1:length(obj$estimated.models))
       {cat("estimation ", i, "\n")
        for (j in 1:length(obj$estimated.models[[i]]$multi.model))
         {r[[i]] <- vector("list",length(obj$estimated.models[[i]]$multi.model))
          cat("model ", j, "\n")
          r[[i]][[j]] <- roots(obj$estimated.models[[i]]$multi.model[[j]])
          if (mod) r[[i]][[j]] <- Mod(r[[i]][[j]])
          print(r[[i]][[j]], digits=digits)
          cat("\n")
         }
       }
    }
  invisible(r)
 }

extract.forecast.cov.estimators.from.model <- function(e,n)
  {# select indicated forecast.cov
   e$forecast.cov <- e$forecast.cov[[n]]
   e$estimation.methods <- e$estimation.methods[[n]]
   e$estimated.models   <- e$estimated.models[[n]]  
   e
}


combine.forecast.cov.estimators.wrt.true <- function(e1,e2)
  {if(! test.equal(e1$true.model, e2$true.model)) 
       warning("true.models are not the same.")
   if(! test.equal(e1$rng == e2$rng)) 
       warning("RNGs are not the same. Second one suppressed.")
   if (!all(c(e1$version[[1]] == e2$version[[1]],
              e1$version[[2]] == e2$version[[2]],
              e1$version[[3]] == e2$version[[3]],
              e1$version[[4]] == e2$version[[4]],
              e1$version[[5]] == e2$version[[5]],
              e1$version[[6]] == e2$version[[6]],
              e1$version[[7]] == e2$version[[7]],
              e1$version[[8]] == e2$version[[8]],
              e1$version[[9]] == e2$version[[9]]))) 
       warning("versions are not the same. Second one suppressed.")
   if(! all(e1$horizon == e2$horizon)) 
       stop("horizon's are not the same.")
   if(e1$discard.before != e2$discard.before) 
       warning("discard.before's are not the same. Second one suppressed.")
   if(e1$est.replications != e2$est.replications) 
       warning("est.replications's are not the same. Second one suppressed.")
   if(e1$pred.replications != e2$pred.replications) 
       warning("pred.replications's are not the same. Second one suppressed.")
   if(! (is.null(e1$simulation.args) & is.null(e2$simulation.args)) )
     if(! all(e1$simulation.args == e2$simulation.args)) 
       warning("simulation.args's are not the same. Second one suppressed.")
   e1$forecast.cov <- append(e1$forecast.cov, e2$forecast.cov)
   e1$estimation.methods <- append(e1$estimation.methods, e2$estimation.methods)
# fix   e1$multi.model <- append(e1$multi.model, e2$multi.model)  
   e1
}

combine.forecast.cov <- function(e1,e2)
  {if(! test.equal(e1$model, e2$model)) 
       warning("models are not the same. Second one suppressed.")
   if(! test.equal(e1$data, e2$data)) 
       warning("data is not the same. Second set suppressed.")
   if(! all(e1$sample.size == e2$sample.size))
       warning("sample.sizes are not the same. Second one suppressed.")
   if(! all(e1$horizons == e2$horizons)) 
       stop("horizon's are not the same.")
   if(e1$discard.before != e2$discard.before) 
       warning("discard.before's are not the same. Second one suppressed.")
   e1$forecast.cov <- append(e1$forecast.cov, e2$forecast.cov)
   e1
}

forecast.cov.reductions.wrt.true <- function(true.model, Spawn=.SPAWN, rng=NULL,
                       simulation.args=NULL,
                       est.replications=2, pred.replications=2,
                       discard.before=10, horizons=1:12,quiet=F,
                       estimation.methods=NULL,
                       criteria=NULL, compiled=.DSECOMPILED)
{# Calculate the forecasts cov of reduced models estimated from simulations of
 # true.model with an estimation method indicated by estimation.methods. 
 #  (estimation.methods is as in estimation.models BUT ONLY THE FIRST IS USED.)
 # discard.before is an integer indicating 1+the number of points in the
 # beginning of forecasts to discard for calculating forecast covariances.
 # criteria can  be a vector of criteria as in information.tests,
 #  (eg c("taic", "tbic") in which case the "best" model for each criteria
 #  is accounted separately. (ie. it is added to the beginning of the list of
 # estimated models)

 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 
 information.criteria <- NULL
 for (i in 1:est.replications)
        {data<-do.call("simulate",append(list(true.model), simulation.args))
         models <-estimate.models(data, trend=T,quiet=quiet,
                       estimation.methods=estimation.methods)
         models$multi.model <-reduced.models.Mittnik(models$multi.model[[1]]) # use only 1
         crit <- NULL
         for (m in models$multi.model) 
            crit<-rbind(crit,information.tests.calculations(l(m, data)))
         if(!is.null(criteria))
           {addmodels <- vector("list", length(criteria))
            for (i in 1:length(criteria))
              addmodels[[i]] <- models$multi.model[[order(crit[,criteria[i]])[1]]]
            models$multi.model <- append(addmodels, models$multi.model)
           }
         information.criteria <- append(information.criteria, list(crit))
         rn <- forecast.cov.wrt.true( models$multi.model, true.model, 
                    pred.replications=pred.replications, zero=T, quiet=quiet,
                    simulation.args=simulation.args, Spawn=Spawn,
                    horizons=horizons, discard.before=discard.before,
                    trend=cbind(1,1:periods(data)) %*% models$trend.coef,
                    compiled=compiled)
         if (i==1)
           r<-rn[c("forecast.cov","forecast.cov.true","forecast.cov.zero","forecast.cov.trend")]
         else
          {for (j in 1:length(models$multi.model))
             r$forecast.cov[[j]] <- r$forecast.cov[[j]]*(i-1)/i + rn$forecast.cov[[j]]/i
           r$forecast.cov.true   <- r$forecast.cov.true*(i-1)/i + rn$forecast.cov.true/i
           r$forecast.cov.zero   <- r$forecast.cov.zero*(i-1)/i + rn$forecast.cov.zero/i
           r$forecast.cov.trend <- r$forecast.cov.trend*(i-1)/i + rn$forecast.cov.trend/i
          }
        }
  classed(append(r, # constructor forecast.cov.estimators.wrt.true (forecast.cov.reductions.wrt.true)
       list(true.model=true.model,
         estimation.methods=c(criteria,estimation.methods),
         rng=rng, version=version,
         horizons=horizons, 
         discard.before=discard.before, est.replications=est.replications,
         pred.replications=pred.replications, simulation.args=simulation.args, 
         information.criteria=information.criteria)),
       c("forecast.cov.estimators.wrt.true", "forecast.cov"))
}


reduced.models.Mittnik <- function(largeModel)
{# Return a list of models with all smaller state dimesions.
  largeModel <- to.SS(largeModel)
  largeModel <- balance.Mittnik(largeModel, n=dim(largeModel$F)[1])
  r <- vector("list", dim(largeModel$F)[1])
  for (j in 1:length(r))
    r[[j]] <- SS(F=largeModel$F[1:j,1:j,drop=F],
                G=if(is.null(largeModel$G)) NULL else largeModel$G[1:j,,drop=F],
                H=largeModel$H[  , 1:j, drop=F],   K= largeModel$K[1:j,,drop=F])
  r
}




############################################################################
#
#       experimental estimation techniques    <<<<<<<<<<
#
############################################################################

est.black.box2 <- function(data, estimation="est.VARX.ls", 
          lag.weight=.9, 
          reduction="reduction.Mittnik", 
          criterion="taic", 
          trend=F, 
          subtract.means=F,  re.add.means=T, 
          standardize=F, verbose=T, max.lag=12)
{if ((estimation!="est.VARX.ls") && (trend) )
     {cat("Trend estimation only support with est.VARX.ls.\n")
      cat("Proceeding using est.VARX.ls.\n")
      estimation<-"est.VARX.ls"
     }

 if(estimation=="est.VARX.ls")
     model <- est.VARX.ls(data,trend=trend, subtract.means=subtract.means, 
                         re.add.means=re.add.means, max.lag=max.lag, 
                         standardize=standardize, lag.weight=lag.weight)
 # else if(estimation=="est.VARX.ar")
 #     model <-est.VARX.ar(data, subtract.means=subtract.means, max.lag=max.lag)
 else
    stop("Only est.VARX.ls estimation technique is supported to date.\n")
 if (verbose) cat("First VAR model,              lags= ", 
          dim(model$model$A)[1]-1, ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(to.SS(model),model$data) # data is standardized if standardize=T in estimation
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model, n= ", 
                  n, ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                     list(model,criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))
    if (verbose) cat("Final reduced state space model, n= ",
              dim(model$model$F)[1], ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose && exists.graphics.device()) check.residuals(model)
 model
}


best.TSestModel     <- function (models, sample.start=10, sample.end=NULL, criterion="aic", verbose=T)
{# return the best model from ... according to criterion
  #  models should be a list of TSestModel's.
  #  models[[i]]$estimates$pred is not recalculated but a sub-sample identified by 
  #  sample.start and  sample.end is used and the likelihood is recalculated. 
  #  If sample.end=NULL data is used to the end of the sample.
  #  taic might be a better default selection criteria but it is not available for ARMA models.
  values <- NULL
  for (lst in models ) 
    {z <- information.tests.calculations(lst, sample.start=sample.start, 
                  sample.end=sample.end)
     values <-rbind(values,z)
    }
  if (verbose)
     {cat("Criterion value for all models based on data starting in period: ",  
           sample.start, "\n")
      cat(values[,criterion], "\n")
     }
#Rbug rbind above looses dimnames
  dimnames(values) <- dimnames(z)
  opt <-order(values[,criterion])[1]  # minimum
  invisible(models[[opt]])
}


est.black.box3 <- function(data, estimation="est.VARX.ls", 
       lag.weight=1.0, 
       reduction="reduction.Mittnik", 
       criterion="aic", 
       trend=F, subtract.means=F,  re.add.means=T, 
       standardize=F, verbose=T, max.lag=12, sample.start=10)
  #  taic might be a better default selection criteria but it is not available for ARMA models.
{if ((estimation!="est.VARX.ls") && (trend) )
     {cat("Trend estimation only support with est.VARX.ls.\n")
      cat("Proceeding using est.VARX.ls.\n")
      estimation<-"est.VARX.ls"
     }
 models <- vector("list", max.lag)
 for (i in 1:max.lag)
   {if(estimation=="est.VARX.ls")
      models[[i]] <- est.VARX.ls(data,trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight)
    else
      stop("Only est.VARX.ls estimation technique is supported to date.\n")
   }
 model <- best.TSestModel(models, criterion=criterion, sample.start=sample.start, verbose=verbose)
 if (verbose) cat("Selected VAR model,              lags= ", 
          dim(model$model$A)[1]-1, ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(to.SS(model),model$data) # data is standardized if standardize=T in estimation
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model,    n= ", 
                  n, ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                      list(model,criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))# , sample.start=sample.start))
    if (verbose) cat("Final reduced state space model, n= ",
              dim(model$model$F)[1], ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose && exists.graphics.device()) check.residuals(model)
 model
}


bft <- function(data, ...) est.black.box4(data, ...)

est.black.box4 <- function(data, estimation="est.VARX.ls", 
                lag.weight=1.0,  variable.weights=1, 
                reduction="reduction.Mittnik", 
                criterion="taic", 
                trend=F, subtract.means=F,  re.add.means=T, 
                standardize=F, verbose=T, max.lag=12, sample.start=10, warn=T)
{if ((estimation!="est.VARX.ls") && (trend) )
     {cat("Trend estimation only support with est.VARX.ls.\n")
      cat("Proceeding using est.VARX.ls.\n")
      estimation<-"est.VARX.ls"
     }
 models <- vector("list", max.lag)
 for (i in 1:max.lag)
   {if(estimation=="est.VARX.ls")
      {model<- est.VARX.ls(data,trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight,
                          warn=warn)
      }
    else if(estimation=="est.wt.variables")
      {model<- est.wt.variables(data, variable.weights,
                        estimation.args=list(trend=trend, 
                          subtract.means=subtract.means,
                          re.add.means=re.add.means, max.lag=i, 
                          standardize=standardize, lag.weight=lag.weight,
                          warn=warn) )
      }
    else
      stop("Estimation technique not yet is supported.\n")
    if (verbose) cat("Estimated  VAR   model       -log likelihood = ", 
        model$estimates$like[1],", lags= ",  dim(model$model$A)[1]-1,"\n")
    model <- l(to.SS(model),model$data,warn=warn) # data is standardized if standardize=T in estimation
    n <- dim(model$model$F)[1]
    if (verbose) cat("Equivalent state space model -log likelihood = ", 
               model$estimates$like[1], ",   n = ", n, "\n")
    if (1 < n)
      {model <- do.call(reduction,
                    list(model,criterion=criterion, verbose=verbose, warn=warn))
      #model <- eval(call(reduction,model,criterion=criterion, 
      #              verbose=verbose, warn=warn))# , sample.start=sample.start))
       if (verbose) cat("Final reduced state space model, n= ",
           dim(model$model$F)[1], ", -log likelihood = ", 
                         model$estimates$like[1], "\n")
      }
     models[[i]] <- model
   }
 model <- best.TSestModel(models, criterion=criterion, sample.start=sample.start, verbose=verbose)
 if (verbose && exists.graphics.device()) check.residuals(model)
 model
}

#  z<-est.black.box4(eg.data,  max.lag=3 ) 



############################################################################
#
#       procedure for testing functions   <<<<<<<<<<
#
############################################################################




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

  invisible(all.ok)
}




dse3.graphics.tests <- function(verbose=T, synopsis=T,  pause=F)
{ if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (synopsis & !verbose) cat("dse3 graphics tests ...")
  if (verbose) cat("  dse3 graphics test 1 ...")
  # If no device is active then write to postscript file 
  if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",
                   width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
            {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
  if(pause) dev.ask(ask=T)

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


############################################################################
#
#       end
#
############################################################################

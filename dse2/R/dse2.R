 
############################################################################

#    functions for data analysis      <<<<<<<<<<<<<

############################################################################


phase.plots <- function(data, max.lag=1,diff=F){
# data is a matrix with a variable in each column.(each column is a
# time series), or an object of class TSdata, inwhich case output.data(data) is used.
# trace plots of data and lagged (difference) data (phase space).
# Non-linearities may show up as a non-linear surface, but this is
#   a projection so, for example, a spherical space would not show up.
#  Some sort of cross-section window would show this but require even more plots.
# A good statistical test would be better!

par(mfcol=c(5,5))  #,mar=c(2.1,4.1,4.1,0.1) )
if (is.TSdata(data)) data <- output.data(data)
Time <- dim(data)[1]
p <-dim(data)[2]
d <- array(NA,c(Time,p,1+max.lag))
d[,,1] <- data
for (l in 1:max.lag) 
  if (diff) {d[(1+l):Time,,l+1] <- d[(1+l):Time,,l]-d[l:(Time-1),,l]}
  else      {d[(1+l):Time,,l+1] <- d[l:(Time-1),,l] } #lag
for (l in 0:max.lag)
   for (l2 in 0:max.lag)
      for(i in 1:p) 
         for(i2 in 1:p) 
           {plot(d[,i,(l+1)],d[,i2,(l2+1)],type="l",xlab="",ylab="")
            title(main=paste("[",i,"L",l,",",i2,"L",l2,"]"))
           }
invisible()
}


############################################################################

#  functions for model estimation (see also VARX in dse.s) and reduction   <<<<<<<<<<<<<

############################################################################


est.wt.variables <- function(data, variable.weights,
                        estimation="est.VARX.ls", estimation.args=NULL)
{ if (is.matrix(variable.weights))
    {if (any(svd(variable.weights)$d == 0))  
       stop("variable.weights transformation must be non singular.")
    }
 else
   {if (any(variable.weights == 0))  stop("variable.weights must be non zero.")
    variable.weights <- diag(variable.weights)
   }
 inv.wts <- solve(variable.weights)
 dimnames(inv.wts)          <-list(NULL, output.series.names(data))
 dimnames(variable.weights) <-list(NULL, output.series.names(data))
 scaled.model <- do.call(estimation, append(list(
           freeze(scale(data, list(output=inv.wts)))), estimation.args))
 model <-scale(TSmodel(scaled.model), list(output=variable.weights))
 model$description <- 
    paste("model estimated by est.wt.variables with", estimation)
 l(model, data)
}

est.min.sqerror <- function(data, init.model, error.weights=1, ...) 
{minimize.TSestModel( l(init.model, data), obj.func=sum.sqerror,
    obj.func.args=list(error.weights=error.weights), ...)
}


est.max.like <- function(emodel, algorithm="nlm",
      max.iter=20, ftol=1e-5, gtol=1e-3,
      dfunc=numerical.grad, line.search="nlm", 
      obj.func=like, obj.func.args=NULL) 
{# maximum likelihood estimation...
 # emodel is an object of class TSestModel (with initial parameter values and data). 
 # max.iter is an integer indicating the maximum number of iterations.
 # The value returned is an object of class TSestModel with  additional
 #  elements $converged, which is T or F indicating convergence, 
 #  and $dfpMin.results or $nlmin.results.
 # If this function calls dfp the Hessian,etc are return as  $dfpMin.results
 # If this function is called again and those results are
 # available then they are used. 
 # This could cause problems if $model is modified. If that is
 # done then $dfpMin.results should be set to NULL.
 # algorithm in {"nlm", "dfpMin", "nlmin", "nlsimplex"}

stop("The current version of est.max.like is superficially broken.")
 if(!is.TSestModel(emodel)) TS.error.exit()
 est  <- emodel$estimates
 Shape <- emodel$model
 Data <- freeze(emodel$data)
 global.assign("Obj.Func.ARGS" , append(list(model=Shape, data=Data),
        obj.func.args))
 if (algorithm=="dfpMin")
    {if (is.null(emodel$dfpMin.results)) parms <- Shape$parms
     else parms <- emodel$dfpMin.results
     dfpMin.results <- dfpMin(obj.func, parms, dfunc=dfunc, 
         max.iter=max.iter, ftol=ftol, gtol=gtol, line.search=line.search) 
     Shape$parms <- dfpMin.results$parms
     emodel <- l(set.arrays(Shape),Data)
     emodel$dfpMin.results <- dfpMin.results
     emodel$converged <- dfpMin.results$converged
     emodel$model$description <- paste("Estimated with max.like/dfpMin (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
    }
 else if (algorithm=="nlmin")
   {nlmin.results <-nlmin(obj.func,Shape$parms, max.iter=max.iter, max.fcal=5*max.iter, ckfc=0.01)
    Shape$parms <- nlmin.results$x
    emodel <- l(set.arrays(Shape),Data)
    emodel$nlmin.results <- nlmin.results
    emodel$converged <- nlmin.results$converged  # this should be improved with conv.type info
    emodel$model$description <- paste("Estimated with max.like/nlmin (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
   }
 else if (algorithm=="nlm")
   {nlm.results <-nlm(obj.func,Shape$parms, hessian=T, iterlim=max.iter)
    Shape$parms <- nlm.results$estimate
    emodel <- l(set.arrays(Shape),Data)
    emodel$nlm.results <- nlm.results
    emodel$converged <- (nlm.results$code<=2)  
    emodel$model$description <- paste("Estimated with max.like/nlm (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
   }
 else if (algorithm=="nlsimplex")
   {results <-nlsimplex(obj.func,Shape$parms, max.iter=max.iter)
    Shape$parms <- results$x
    emodel <- l(set.arrays(Shape),Data)
    emodel$nlsimplex.results <- results
    emodel$converged <- results$converged  # this should be improved with conv.type info
    emodel$model$description <- paste("Estimated with max.like/nlsimplex (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
   }
  else stop(paste("Minimization method ", algorithm, " not supported."))
# remove(c("Obj.Func.ARGS"), where=1)
 emodel
}    


est.black.box <- function(data,...)
{# call current best black box technique.
  est.black.box4(data, ...)
}


est.black.box1 <- function(data,estimation="est.VARX.ls", reduction="reduction.Mittnik", 
        criterion="taic", trend=F, subtract.means=F, verbose=T, max.lag=6)
{if ((estimation!="est.VARX.ls") && (trend) )
     {cat("Trend estimation only support with est.VARX.ls.\n")
      cat("Proceeding using est.VARX.ls.\n")
      estimation<-"est.VARX.ls"
     }

 if(estimation=="est.VARX.ls")
     model <- est.VARX.ls(data,trend=trend, subtract.means=subtract.means, 
                          max.lag=max.lag)
 else if(estimation=="est.VARX.ar")
     model <- est.VARX.ar(data, subtract.means=subtract.means, max.lag=max.lag)
 else if(estimation=="est.VARX.ls")
     model <- est.VARX.ls(data,trend=trend, subtract.means=subtract.means, 
                        max.lag=max.lag)
 else if(estimation=="est.SS.Mittnik")
     model <- est.SS.Mittnik(data,max.lag=max.lag, 
                             subtract.means=subtract.means, normalize=F)
 else
   stop("estimation technique not supported.")
 if (verbose) 
   cat("First VAR model,              lags= ", dim(model$model$A)[1]-1,
       ", -log likelihood = ", model$estimates$like[1], "\n")
 model <- l(to.SS(model),data)
 n <- dim(model$model$F)[1]
 if (verbose) cat("Equivalent    state space model, n= ", n,
                  ", -log likelihood = ", model$estimates$like[1], "\n")
 if (1 < n)
   {model <- do.call(reduction,
                     list(model, criterion=criterion, verbose=verbose))
   #model <- eval(call(reduction,model,criterion=criterion, verbose=verbose))
    if (verbose) 
       cat("Final reduced state space model, n= ", dim(model$model$F)[1],
           ", -log likelihood = ", model$estimates$like[1], "\n")
   }
  if (verbose && exists.graphics.device()) check.residuals(model)
 model
}


est.SS.Mittnik <- function(data, max.lag=6, n=NULL, subtract.means=F, normalize=F)
{#  estimate a nested-balanced state space model by svd from least squares
 # estimate of markov parameters a la Mittnik p1195 Comp.Math Appl.v17,1989.
 # The quality of the estimate seems to be quite sensitive to max.lag, 
 #   and this is not properly resolved yet.
 # If n is not supplied the svd criteria will be printed and n prompted for.
 # If subtract.means=T then the sample mean is subtracted. 
 # If normalize is T the lsfit estimation is done with outputs normalize to cov=I
 # (There still seems to be something wrong here!!).
 # The model is then re-transformed to the original scale.
 
  data <- freeze(data)
  m <- ncol(input.data(data))
  if(is.null(m))  m <- 0
  p <- ncol(output.data(data))
  N <- nrow(output.data(data))
  if (subtract.means)
    {if(m!=0)input.data(data)<-input.data(data)-t(matrix(apply(input.data(data),2, mean), m,N))
     output.data(data)<- output.data(data) - t(matrix(apply(output.data(data),2, mean), p,N))
    }
  if (normalize)
    {svd.cov <- svd(var(output.data(data)))
     output.data(data) <- output.data(data) %*% svd.cov$u %*% diag(1/svd.cov$d^.5)
    }
      # shift input to give feedthrough in one period
  if (m != 0) {z <- cbind(input.data(data)[2:N,],output.data(data)[1:(N-1),])}
  else z <- output.data(data)
  Past <- matrix(NA,N-1-max.lag,(p+m)*(1+max.lag))
  for (i in 0:max.lag) 
    Past[,(1+(m+p)*i):((m+p)*(1+i))] <-z[(1+max.lag-i):(N-1-i),]
  M <- t(lsfit(Past,output.data(data)[(max.lag+2):N,],intercept=F)$coef)
  if (normalize && (m!=0))  # correct exogenous blocks for normalization
    {Tinv <- diag(svd.cov$d^.5)%*%svd.cov$u
     for (i in 0:max.lag) 
       M[,(1+(m+p)*i):(m+(m+p)*i)] <- Tinv %*% M[,(1+(m+p)*i):(m+(m+p)*i),drop=F]
    }
  if (p==1) M <-matrix(M,1,length(M))
#  browser()
  model <-balance.Mittnik.svd( M, m=m, n=n )$model
  z <-"nested-balanced model from least sq. estimates of markov parameters a la Mittnik"
  if(subtract.means) z <-paste(z," - means subtracted")
  if(normalize)      z <-paste(z," - outputs normalized")
  model$description  <-z
  series.names(model)  <- series.names(data)
  l(model, data)
}


reduction.Mittnik <- function(model, data=NULL, criterion=NULL, verbose=T,warn=T)
{# nested-balanced state space model reduction by svd of Hankel generated from a model
# If a state space model is supplied the max. state dimension for the result is
#  taken from the model. If an ARMA model is supplied then singular values 
#  will be printed and the program prompts for the max. state dimension.

  if(!is.TSm.or.em(model)) TS.error.exit()
  if (is.TSestModel(model)) 
    {if (is.null(data)) data <-TSdata(model)
     model <- TSmodel(model)
    }
  if(is.null(data))
    stop("Reduction requires data to calculate criteria (balancing does not).")
  nMax <- if(is.SS(model)) dim(model$F)[2] else NULL
  reduction.Mittnik.from.Hankel( markov.parms(model), nMax=nMax, data=data, 
        criterion=criterion, verbose=verbose, warn=warn)
}

reduction.Mittnik.from.Hankel <- function(M, data=NULL, nMax=NULL, criterion=NULL, verbose=T, warn=T, spawn=.SPAWN)
{# Select a reduced state space model by svd a la Mittnik.
 #  Models and several criteria for all state dimensions up to the max.  
 #  state dim. specified are calculated. (If nMax is not supplied then svd
 #  criteria are printed and the program prompts for nMax). 
 # The output dimension p is taken from nrow(M).
   # M is a matrix with p x (m+p)  blocks giving the markov parameters,
  # that is, the first row of the Hankel matrix. It can be generated from the
  # model as in the function markov.parms, or from the data, as in the function
  # est.SS.Mittnik.

  # data is necessary only if criteria (AIC,etc) are to be calculated.

 # See the documentation for reduction.Mittnik.

   data <- freeze(data)
   m <-ncol(input.data(data))      # dim of input series
   if(is.null(m))m<-0
   z <-balance.Mittnik.svd(M, m, nMax)
   largeModel <- z$model
   svd.crit    <-z$crit
   n <- dim(largeModel$F)[1]
   if (!spawn)
     {# The more complicated For loop used below is to avoid S memory problems
      #  with the more straight forward version:
      values <- NULL 
      for (i in 1:n) 
        {if(m!=0) z <-largeModel$G[1:i,,drop=F]
         else     z <-NULL
         z <-SS(F=largeModel$F[1:i,1:i,drop=F],G=z,
                  H=largeModel$H[,1:i,drop=F],K= largeModel$K[1:i,,drop=F])
         z <-information.tests.calculations(l(z,data, warn=warn))
         values <-rbind(values, z)
         if (verbose) cat(".")
        }
      }
   else 
     {if (verbose) cat("Spawning processes to calculate criteria test for state dimension 1 to ",n)
      forloop <- function(largeModel, data, warn=T)
           {if(!is.null(largeModel$G)) z <-largeModel$G[1:forloop.i,,drop=F]
            else                       z <-NULL
            z <-SS(F=largeModel$F[1:forloop.i,1:forloop.i,drop=F],G=z,
                     H=largeModel$H[  , 1:forloop.i, drop=F],
                     K=largeModel$K[1:forloop.i,,drop=F])
            information.tests.calculations(l(z,data, warn=warn))
           }
       assign("balance.forloop", forloop, where=1)
       assign("forloop.n", n,   where=1 )
       assign("forloop.values", matrix(NA,n,12),   where=1 )
       assign("forloop.largeModel", largeModel, where=1)
       assign("forloop.data", data,   where=1 )
       assign("forloop.warn", warn,   where=1 )
       on.exit(remove(c("balance.forloop", "forloop.i", "forloop.n", 
          "forloop.values", "forloop.largeModel", 
          "forloop.data", "forloop.warn"),where=1))
       For (forloop.i=1:forloop.n,
         forloop.values[forloop.i,]<-balance.forloop(forloop.largeModel, 
         forloop.data, forloop.warn), sync=T)
       values <-forloop.values
      }
    dimnames(values) <- list(NULL,c("port","like","aic","bic", 
          "gvc","rice","fpe","taic","tbic","tgvc","trice","tfpe")) 
    if (verbose) cat("\n")
    opt <-apply(values,2,order)[1,]  # minimum
    if (verbose | is.null(criterion))
      {zz<-criteria.table.nheading()
       options(width=120)
       print(values,digits=4)
       cat("opt     ")
       for (i in 1:length(opt)) cat(opt[i],"    ")
       cat("\n")
       zz<-criteria.table.legend()
      }
    if (is.null(criterion))
      {n <- read.int("Enter the state dimension (enter 0 to stop): ")
       if( n<1) stop("TERMINATED! STATE DIMENSION MUST BE GREATER THAN 0!")
      }
    else { n <- opt[criterion == dimnames(values)[[2]]]  }
    if(m==0) z <-NULL 
      else   z <-largeModel$G[1:n,,drop=F]
    model <- SS(description="nested model a la Mittnik",
          F=largeModel$F[1:n,1:n,drop=F],G=z,
          H=largeModel$H[,1:n,drop=F],K= largeModel$K[1:n,,drop=F], 
          names=series.names(data))
    l(model,data, warn=warn)          
}


############################################################################

#    functions for model analysis   <<<<<<<<<<<<<

############################################################################


shock.decomposition <- function(model, horizon=30, shock=rep(1,horizon))
{ if(!is.TSm.or.em(model)) TS.error.exit()
  m <-input.dimension(model)
  p <-output.dimension(model)

  if (is.TSestModel(model)) 
     {data <- TSdata(model)   # just for input
      model <- TSmodel(model)
     }
  else 
    {if( m > 0 ) data <- TSdata(input=matrix(0,horizon,m))
     else        data <- TSdata(output=matrix(0,horizon,p))
    }
   model$z0 <- NULL    # zero initial conditions
   par(mfrow = c(p, p) , mar = c(2.1, 4.1,3.1, 0.1) )
   for (i in 1:p)
     {output.data(data) <- matrix(0, horizon, p)
      output.data(data)[,i] <- shock   
      z <- l(model,data)
      tfplot(z, reset.screen=F)
     }
invisible()
}

############################################################################

#    functions for forecasting    <<<<<<<<<<<<<

# Class "feather.forecasts" has a forecast path from multiple starting points
#  in the data (so the graph may look like a feather).
# In the simplest case it would start from the end of the data 
#  and give the path out to a horizon.

############################################################################
############################################################################

#    methods for forecast       <<<<<<<<<<<<<

############################################################################

is.forecast <- function(obj) inherits(obj,"forecast")

forecast <- function(obj, ...)   UseMethod("forecast")
forecast.TSestModel <- function(obj, ...){forecast(TSmodel(obj),TSdata(obj),...)}
forecast.TSdata <- function(obj, model, ...){forecast(model, obj, ...)}

forecast.TSmodel <- function(obj, data,  horizon=36, conditioning.inputs=NULL, conditioning.inputs.forecasts=NULL, percent=NULL)
{# obj must be a TSmodel
 # Calculate (multiple) forecasts from the end of data to a horizon determined
 # either from supplied input data or the argument horizon (more details below).
 # In  the case of a model with no inputs the horizon is determined by the
 #   argument horizon.
 # In the case of models with inputs, on which the forecasts
 #  are conditioned, the argument horizon is ignored (except when percent is
 #  specified) and the actual horizon is determined by the inputs in the 
 #  following way:
 # If inputs are not specified by optional arguments (as below) then the default
 #  will be to use input.data(data). This will be the same as the function l() unless
 #  input.data(data) is longer (after NAs are trimmed from each separately) than
 #  output.data(data).
 # Otherwise, if conditioning.inputs is specified it is used for input.data(data).
 #    It must be a time series matrix or a list of time series matrices each
 #    of which is used in turn as input.data(data). The default above is the same as
 #        forecast(model, trim.na(data), conditioning.inputs=trim.na(input.data(data)) )
 # Otherwise, if conditioning.inputs.forecasts is specified it is appended 
 #   to input.data(data). It must be a time series  
 #   matrix or a list of time series matrices each of which is 
 #   appended to input.data(data) and the concatenation used as conditioning.inputs.
 #   Both conditioning.inputs and conditioning.inputs.forecasts should not be
 #   specified.
 # Otherwise, if percent is specified then conditioning.inputs.forecasts are set
 #    to percent/100 times the value of input corresponding to the last period
 #    of output.data(data) and used for horizon periods. percent can be a vector, 
 #    in which case each value is applied in turn. ie c(90,100,110) would would 
 #    give results for conditioning.input.forecasts 10 percent above and below 
 #    the last value of input.

 # The result is an object of class forecast which is a list with 
 #   elements model, horizon, conditioning.inputs, percent, and forecast.
 #   forecast is a list with TSdata objects as elements, one for each element 
 #   in the list conditioning.inputs.
 
 if(!is.TSmodel(obj)) stop("obj must be a TSmodel in forecast.TSmodel")

 output <- trim.na(output.data(data))
 sampleT <- dim(output)[1]

 if (0==input.dimension(obj))
     {if (0 != (input.dimension(data)))
              warning("data has input and model does not take inputs.")
      if (!is.null(conditioning.inputs))
          warning("model does not take inputs. conditioning.inputs ignored.")
      if (!is.null(conditioning.inputs.forecasts))
   warning("model does not take inputs. conditioning.inputs.forecasts ignored.")
      if (!is.null(percent))
          warning("model does not take inputs. percent ignored.")
      pr <- l(obj, data, sampleT =sampleT, 
                           predictT=sampleT+horizon)$estimates$pred
      pred <- tfwindow(pr, end=end(output), warn=F)
      pr <- tfwindow(pr, start=c(0,1)+end(output), warn=F)
    #  pr[1:(sampleT-1),] <- NA
    #  pr[sampleT,] <- output[sampleT,]
      proj <- list(pr)
     }
 else
  {if ((!is.null(conditioning.inputs)) &
       (!is.null(conditioning.inputs.forecasts)))
       warning(paste("conditioning.inputs and conditioning.inputs.forecasts",
        " should not both be supplied. conditioning.inputs are being used."))

   if ((!is.null(conditioning.inputs))& (!is.null(percent)))
       warning(paste("conditioning.inputs and percent",
         " should not both be supplied. conditioning.inputs are being used."))

   if ((!is.null(percent))& (!is.null(conditioning.inputs.forecasts)))
      warning(paste("percent and conditioning.inputs.forecasts should not",
          " both be supplied. conditioning.inputs.forecasts are being used."))

   if (!is.null(conditioning.inputs)) {} # do nothing
   else if (!is.null(conditioning.inputs.forecasts))
         {if (is.matrix(conditioning.inputs.forecasts)) 
            conditioning.inputs.forecasts <-list(conditioning.inputs.forecasts)
          conditioning.inputs <- list()
          for (i in 1:length(conditioning.inputs.forecasts) )
            {inp <-tframed(rbind(input.data(data),conditioning.inputs.forecasts[[i]]), 
              list(start=start(input.data(data)),
                   frequency=frequency(input.data(data))))
             conditioning.inputs <- append(conditioning.inputs, list(inp))
         }  }  
   else if (!is.null(percent))   
        {last.in <- input.data(data)[sampleT,]
         for (i in 1:length(percent) )
           {pol <- t(matrix(last.in*percent[i]/100, length(last.in), horizon))
            inp <-ts(rbind(input.data(data)[seq(sampleT),,drop=F],pol), 
                     start=start(input.data(data)),
                     frequency=frequency(input.data(data)))
            conditioning.inputs <- append(conditioning.inputs, list(inp))
        }  }  
   else conditioning.inputs <- trim.na(input.data(data))

   if (is.matrix(conditioning.inputs))
          conditioning.inputs <- list(conditioning.inputs)

   proj <- NULL
   for (policy in  conditioning.inputs)
        {pdata <- TSdata(input=policy, output=output)
         if(2 != length(start(pdata)))
            stop("input and output data must have the same starting period (after NAs are removed).")
      predictT <- dim(policy)[1]
      if (sampleT > predictT) 
         stop("input series must be at least as long as output series (after NAs are removed).")
      horizon <- predictT - sampleT
      pr <- l(obj, pdata, sampleT=sampleT, predictT=predictT)$estimates$pred
#    The following lines sometimes cause problems if output is output[...,]
#    See comments in dse2.function.tests
        if (0 == length(proj)) pred <- tfwindow(pr, end=end(output), warn=F)
        if(all(end(pr)==end(output)))
          {pr <- NULL
           warning("Input is not longer than output data. No forecasts produced.") 
          }
        else pr <- tfwindow(pr, start=c(0,1)+end(output), warn=F)
     #    pr[1:(sampleT-1),] <- NA
     #    pr[sampleT,] <- output[sampleT,]  # so plots show first step
         proj <- append(proj, list(pr))
        }
   }
 invisible(classed(list(model=obj, data=data,  # forecast constructor
                horizon=horizon, percent=percent,
                conditioning.inputs=conditioning.inputs,
                conditioning.inputs.forecasts=conditioning.inputs.forecasts,
                forecast=proj,pred=pred), "forecast"))
}


# extract the forecasts
forecasts <- function(obj, ...)   UseMethod("forecasts")
forecasts.forecast <- function(obj, ...){obj$forecast}


test.equal.forecast <- function(obj1, obj2, fuzz=1e-14)
{# N.B. models are not compared (so equivalent models will compare true)
 # inputs are not compared, so they may be provided differently.
 r <- all(dseclass(obj1) == dseclass(obj2))
 if (r) r <- all(output.data(obj1$data) == output.data(obj2$data))
 if (r) r <- all(obj1$horizon == obj2$horizon)
 if (r) r <- fuzz > max(abs(obj1$pred - obj2$pred))
 if (r) r <- length(obj1$forecast)==length(obj2$forecast)
 for (i in seq(length(obj1$forecast)))
   if (r) r <- fuzz > max(abs(obj1$forecast[[i]] - obj2$forecast[[i]]))
 r
}

tfplot.forecast <- function(x, start.=NULL, end.=NULL,
        select.series = seq(length=output.dimension(x$data)),
        names = output.series.names(x$data))
{#The default starting point (start.) for plots is the start data.
 #The default ending point (end.) for plots is the end of forecast.
   if (is.null(x$forecast[[1]]))
      stop("Object to be plotted contains no forecast")
   output <-trim.na(output.data(x$data))
   old.par <-par(mfcol = c(length(select.series), 1), mar= c(5.1,6.1,4.1,2.1))
   on.exit(par(old.par))
   N <- length(x$forecast)
   H <- 0
   for (t in 1:N) H <- max(H, dim(x$forecast[[t]])[1])
   tf <-expand(tframe(output), add.end=H)
   for(i in select.series) 
        {z <- c(output[,i], rep(NA,H))
         for (t in 1:N)
            {zz <- c(rep(NA,periods(output)),x$forecast[[t]][,i],
                     rep(NA,H-dim(x$forecast[[t]])[1]))
             zz[periods(output) ] <- output[periods(output), i] #so line joins last data to first forecast
             z <- cbind(z,zz)
            }
         tframe(z) <- tf
         if (!is.null(start.)) z <- tfwindow(z,start=start., warn=F)
         if (!is.null(end.))   z <- tfwindow(z,end=end., warn=F)
         tfplot(z, ylab = names[i])
         if(i == select.series[1]) 
             title(main = "Predictions (dotted) and actual data (solid)")
        }
   invisible()
}

output.series.names.forecast <- function(obj)
   {m <- output.series.names(obj$model)
    d <- output.series.names(obj$data)
    if(!all(m == d))
       warning("data and model names do not correspond. Model names returned.")
    m
   }

input.series.names.forecast <- function(obj)
   {m <- input.series.names(obj$model)
    d <- input.series.names(obj$data)
    if(!all(m == d))
       warning("data and model names do not correspond. Model names returned.")
    m
   }

############################################################################

#    methods for feather.forecasts        <<<<<<<<<<<<<

############################################################################

is.feather.forecasts <- function(obj) inherits(obj, "feather.forecasts")

output.series.names.feather.forecasts <- function(x) output.series.names(x$data)
 input.series.names.feather.forecasts <- function(x)  input.series.names(x$data)

feather.forecasts <- function(obj, ...) UseMethod("feather.forecasts")

feather.forecasts.TSestModel <- function(obj, data=NULL, ...)
     {if (is.null(data)) data <- TSdata(obj)
      feather.forecasts(TSmodel(obj), data, ...)}

feather.forecasts.TSdata <- function(obj, model, ...)
     {feather.forecasts(model, obj, ...)}

feather.forecasts.TSmodel <- function(model, data, horizon=36,
             from.periods =NULL, ...)
  {if(!is.TSmodel(model)) TS.error.exit(clss="TSmodel")
   if(!is.TSdata(data)) TS.error.exit(clss="TSdata")
   data <- freeze(data)
   if (is.null(from.periods))
     {if(0 == output.dimension(data)) from.periods <-
             10*seq(floor(periods(data)/10))
      else from.periods <-
             10*seq(floor(min(periods(data), input.periods(data)-horizon)/10))
     }
   # periods.TSPADIdata returns NA rather than fetching data. Note:Previously freeze
   #  was not done above and pr below just left NULL for TSPADIdata, so some
   #  of this can be cleaned out.
   if ((!is.na(periods(data))) && (max(from.periods) > periods(data) ))
     stop("from.periods cannot exceed available output data.")
   if (0 != (input.dimension(data)))
     if ((!is.na(input.periods(data))) && 
        ((max(from.periods)+horizon) > input.periods(data) ))
       stop("forecasts cannot exceed available input data.")
   shf <- start.shift(model,data,y0=NULL)  # ? y0=y0)
   proj <- NULL
   for (sampleT in from.periods)
     {pr <-l(model, data, sampleT=sampleT, 
              predictT=sampleT+horizon, result="pred", ...)
      pr[1:(sampleT-1),] <- NA
      # make period before prediction = data so graphics look correct.
      #following if is kludge to prevent retrieving data
      pr[sampleT,] <- output.data(data)[sampleT+shf$shift*shf$lags,]
      proj <- append(proj, list(pr))
     }
   # names are available from  data or model
   invisible(classed(list(model=model, # feather.forecasts constructor
                data=data, from.periods=from.periods, 
                horizon=horizon, feather.forecasts=proj), "feather.forecasts"))
}


forecasts.feather.forecasts <- function(obj){obj$feather.forecasts}

tfplot.feather.forecasts <- function(x, start.=NULL, end.=NULL, select.series=NULL, graphs.per.page=5, reset.screen=T)
{#The default starting point (start.) for plots is the start of data.
 #The default ending point (end.) for plots is the end of forecasts.
   p <- dim(x$feather.forecasts[[1]])[2]
   freq <- frequency(x$data)
   names <- output.series.names(x)
   if(is.null(names)) names <- paste("output", 1:p)
   if (is.null(start.)) start. <- start(x$data)
   if (is.null(end.))   end.   <- add.date(end(output.data(x$data)),
                                   max(x$horizon), frequency(x$data))
   if (is.null(select.series)) select.series <- 1:p
   if (!is.numeric(select.series)) select.series <- match(select.series, names)
   if(reset.screen) 
     {Ngraphs <- length(select.series)
      Ngraphs <- min(Ngraphs, graphs.per.page)
      old.par <- par(mfcol = c(Ngraphs, 1), mar= c(5.1,6.1,4.1,2.1)) 
      on.exit(par(old.par))
     }
   # if below is a kludge to skip getting TSPADI data.  
   for(i in select.series) 
        {#if (is.TSPADIdata(x$data)) 
         #  {zz <- NULL # kludge
         #   ltys <- rep(2,length(x$from.periods))
         #  }
         #else 
           {zz <- tfwindow(output.data(x$data,series=i), start=start.,warn=F)
            ltys <- c(1,rep(2,length(x$from.periods)))
           }
         for (t in 1:length(x$from.periods))
            {zz <- tbind(zz,
                     ts(x$feather.forecasts[[t]][,i], 
                       start=start(x$feather.forecasts[[t]]),freq=freq))
            }
         tfplot(tfwindow(zz,start=start.,end=end., warn=F), ylab=names[i], lty=ltys)
         if(i == select.series[1]) 
             title(main = "Predictions (dotted) and actual data (solid)")
        }
   invisible()
}



############################################################################
#
#       procedure for testing functions   <<<<<<<<<<<<<
#
############################################################################



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
  z <- feather.forecasts( modSS,  from.periods=c(250,300))
  error <- max(abs
       (c(z$feather.forecasts[[1]][286,],z$feather.forecasts[[2]][336,])
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

  invisible(all.ok)
}

dse2.graphics.tests <- function(verbose=T, synopsis=T,  pause=F)
{# graphics tests do not do any value comparisons
  if (synopsis & !verbose) cat("dse2 graphics tests ...")
  
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (verbose) cat("  dse2 graphics test 1 ...")

  # If no device is active then write to postscript file 
  if (!exists.graphics.device())
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=F, print.it=F, append=F)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }
  if(pause) dev.ask(ask=T)

  data <- eg1.DSE.data.diff
  mod1 <- TSmodel(est.VARX.ls(data,max.lag=3))
  modSS <- l(to.SS(mod1),data)

  z <- feather.forecasts( modSS,  from.periods=c(230,250))
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


############################################################################
#
#       end
#
############################################################################

 
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

#  functions for model estimation (see also VARX ) and reduction   <<<<<<<<<<<<<

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


est.max.like <- function(Shape, ...) {UseMethod("est.max.like")}

est.max.like.TSestModel <- function(Shape, data=TSdata(Shape), ...) {
	 # if Shape is result from a previous est.max.like then the gradient
	 # hessian and other information should be extracted, but
	 est.max.like(TSmodel(Shape), data, ...) }

est.max.like.TSdata <- function(data, Shape, ...) {
	 est.max.like(data, TSmodel(Shape), ...) }

est.max.like.TSmodel <- function(Shape, data, 
	algorithm="optim",
	algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf, hessian=TRUE)
	)
{# maximum likelihood estimation
 # "nml" algorithm.args=list(hessian=T, iterlim=20, 
 #     dfunc=numerical.grad, line.search="nlm",ftol=1e-5, gtol=1e-3,)
 data <- freeze(data)
 func.like <- function(parms,Shape,data)
      {l(set.arrays(Shape,parms=parms),data,result="like") }

 if (algorithm=="optim")
    {results <- optim(Shape$parms, func.like, method=algorithm.args$method,
	gr=algorithm.args$gr, 
	lower=algorithm.args$lower, upper=algorithm.args$upper,
	control=algorithm.args$control, hessian=algorithm.args$hessian,
	Shape, data) 
     emodel <- l(set.arrays(Shape, parms=results$par),data)
     emodel$est$algorithm <- algorithm
     emodel$est$results <- results
     emodel$est$converged <- results$converged
     emodel$model$description <- paste("Estimated with max.like/optim (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
    }
 else if (algorithm=="nlm")
   {warning("This has not been tested recently (and there have been changes which may affect it.")
    results <-nlm(func.like,Shape$parms, hessian=algorithm.args$hessian, 
    	iterlim=algorithm.args$iterlim)
    emodel <- l(set.arrays(Shape, parms=results$estimate),data)
    emodel$est$algorithm <- algorithm
    emodel$est$results <- results
    emodel$est$converged <- results$code <= 2
    emodel$model$description <- paste("Estimated with max.like/nlm (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
   }
 else if (algorithm=="nlmin")
   {warning("This has not been tested recently (and there have been changes which may affect it.")
     results <-nlmin(func.like,Shape$parms, max.iter=algorithm.args$max.iter, 
     	max.fcal=5*algorithm.args$max.iter, ckfc=0.01)
     emodel <- l(set.arrays(Shape, parms=results$parms),data)
     emodel$est$algorithm <- algorithm
     emodel$est$results <- results
     emodel$est$converged <- results$converged
    # above should be improved with conv.type info
    emodel$model$description <- paste("Estimated with max.like/nlmin (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
   }
 else if (algorithm=="dfpMin")
    {stop("This optimization method is no longer supported.")
     results <- dfpMin(func.like, Shape$parms, 
	dfunc=algorithm.args$dfunc, 
	max.iter=algorithm.args$max.iter, 
	ftol=algorithm.args$ftol, gtol=algorithm.args$gtol, 
	line.search=algorithm.args$line.search) 
     emodel <- l(set.arrays(Shape, parms=results$parms),data)
     emodel$est$algorithm <- algorithm
     emodel$est$results <- results
     emodel$est$converged <- results$converged
     emodel$model$description <- paste("Estimated with max.like/dfpMin (",
       c("not converged", "converged")[1+emodel$converged],
       ") from initial model: ", emodel$model$description)
    }
  else stop(paste("Minimization method ", algorithm, " not supported."))
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
  if (verbose &&  dev.cur() != 1 ) check.residuals(model)
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
         tfOnePlot(z, ylab = names[i])
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
         tfOnePlot(tfwindow(zz,start=start.,end=end., warn=F), ylab=names[i], lty=ltys)
         if(i == select.series[1]) 
             title(main = "Predictions (dotted) and actual data (solid)")
        }
   invisible()
}




############################################################################

# Functions in the next group are mainly for evaluating estimation techniques.

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

monte.carlo.simulations.default <- function (model, simulation.args = NULL, 
 		replications = 100, rng = NULL, quiet = F){
	if (is.null(rng)) rng <- set.RNG()
	else {
		old.rng <- set.RNG(rng)
		on.exit(set.RNG(old.rng))
	}
	arglist <- append(list(model), simulation.args)
	r <- do.call("simulate", arglist)
	if (! is.matrix(r)) stop("simulate(model) must return a matrix.")
        result <- array(NA, c(dim(r), replications))
        tfr <- tframe(r)
        result[, , 1] <- r
        if (1 < replications) 
            for (i in 2:replications)
	      result[, , i] <- do.call("simulate", arglist)
	series.names(result) <- series.names(r)
	result <- tframed(result, tfr)
	invisible(classed(list(simulations = result, model = model, 
        rng = rng,  simulation.args = simulation.args, 
        description = "data generated by monte.carlo.simulation.default"), 
        c("monte.carlo.simulations")))
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
        result[,,i] <- output.data(do.call("simulate", arglist))
  }
series.names(result) <- output.series.names(model)
result <- tframed(result, tfr)  # my more general multidimensional ts
invisible(classed( # constructor monte.carlo.simulations
         list(simulations=result, model=model, rng=rng, simulation.args=simulation.args,
              description = "data generated by monte.carlo.simulations.TSmodel"),
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
  {names <- series.names(x$simulations)
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
         tfOnePlot(tfwindow(zz,start=start.,end=end., warn=F), ylab=names[i]) #tsplot
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
        stop("Neither ksmooth nor density are available to calculate the plot.")
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
        stop("Neither ksmooth nor density are available to calculate the plot.")
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


tfplot.parms.ee <- function(x, cum=T, norm=F, bounds=T, invert=F, Sort=F,
	graphs.per.page = 5)
{# if cum is true the cummulative average is plotted.
 # if norm is true the norm is used, each parameter is plotted.
 # if invert is true the reciprical is used (before cummulating).
 # if Sort is true then sort is applied (before cum). This is not usually
 #   recommended but of interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
 # Plotting the true lines only makes sense if truth is the same length as 
 #  result (and sometimes not even then). 
      N<-length(x$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((x$result)[[i]]))
      r <- matrix(0,N,n)
      plottrue <- T
      for (i in 1:N) {
	ni <- length((x$result)[[i]])
	r[i,1:ni] <- (x$result)[[i]]
	if (ni != n) plottrue <- F
	}
      if (invert) r <- 1/r
      if(norm)    r <- matrix((apply(r^2,1,sum))^.5, N,1)
      r[is.infinite(r)] <- 0
      if (Sort) r <- t(apply(r,1,sort))
      if (n != length(x$truth)) plottrue <- F
      if (plottrue) {
	true.lines <- c(x$truth)
	if (invert) true.lines <- 1/true.lines
	if(norm) true.lines <-sum(true.lines^2)^.5
	if (Sort) true.lines <- sort(true.lines)
	true.lines <-t(matrix(true.lines, length(true.lines),N))
	if (bounds) {
		z  <- r-true.lines
		Om <- t(z) %*% z/(nrow(z)-1)
		Om <- diag(Om)^.5
		Om <- t(matrix(Om, length(Om), N))
		Om <- Om/matrix((1:N)^.5 , N, ncol(Om))
        }
      }
      if (cum) r<- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      series.names(r) <- paste("parm", seq(ncol(r)))
#      matplot(x=matrix(seq(nrow(r)),nrow(r),1), y=cbind(0,true.lines,r, Om), 
#              type="l", lty=c(1,rep(3,dim(true.lines)[2]), rep(4,dim(r)[2]), 
#                     rep(2,2*dim(r)[2])))
      if (plottrue & bounds)
         tfplot(r, true.lines,true.lines + Om, true.lines - Om,
            graphs.per.page = graphs.per.page)
      else if (plottrue) tfplot(r, true.lines, graphs.per.page=graphs.per.page)
      else               tfplot(r, graphs.per.page = graphs.per.page)
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

distribution.parms.ee <- function(obj,  Sort=F, bandwidth=0.2,
	graphs.per.page=5)
{# if Sort is true then sort is applied (before cum). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
      N<-length(obj$result)
      n <- length(obj$truth)
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      if (n == length(obj$truth)) plottrue <- T
      else {
         warning("Number of true parameters does not match number estimated.")
	 plottrue <- F
      }
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,1:length((obj$result)[[i]])] <- (obj$result)[[i]]
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      xlab <-"parameter "
      old.par <- par(mfcol=c(min(graphs.per.page, ncol(r)),1))
      par(mar=c(5.1, 4.1, 4.1, 2.1))
      on.exit(par(old.par))
      for ( i in 1:ncol(r))
         {if (is.Splus()) rd <- ksmooth(r[,i], bandwidth=bandwidth)
          if (is.R())     rd <- density(r[,i], bw=bandwidth)
          plot(rd, type="l", ylim=c(0, max(rd$y)),
               ylab="density",  xlab=paste(xlab, i) , main="")
          if (plottrue) lines(rep(true.lines[i],2),c(1,0))
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
                  as.double(rep(0,is)),         # scratch array
                  DUP=TRUE)$proj
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
                  as.integer(gain),
                  as.double(z),
                  as.double(P), 
		  DUP=TRUE)$proj
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
                  as.double(err), 
		  DUP=TRUE) [c("forecast.cov","sample.size")]
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

forecast.cov.TSestModel <- function(obj, data=NULL, ..., discard.before=NULL,
       horizons=1:12, zero=F, trend=F, estimation.sample= NULL,
       compiled=.DSECOMPILED)
 {forecast.cov(TSmodel(obj), ..., data=if(is.null(data)) TSdata(obj) else data,
	discard.before=discard.before,
	horizons=horizons, zero=zero, trend=trend,
	estimation.sample= if(is.null(estimation.sample)) periods(data) else estimation.sample,
	compiled=compiled)}

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
                  as.double(rep(0,is)),         # scratch array
                  DUP=TRUE)[c("forecast.cov","sample.size")]
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
     .Fortran("kfepr",
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
                  as.integer(gain),
                  as.double(z),
                  as.double(P), 
		  DUP=TRUE) [c("forecast.cov","sample.size")]
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
  if (verbose &&  dev.cur() != 1 ) check.residuals(model)
 model
}


best.TSestModel <- function(models, sample.start=10, sample.end=NULL, criterion="aic", verbose=T)
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
  if (verbose &&  dev.cur() != 1 ) check.residuals(model)
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
 if (verbose &&  dev.cur() != 1 ) check.residuals(model)
 model
}

#  z<-est.black.box4(eg.data,  max.lag=3 ) 



# Functions in the next group are mainly for evaluating the information <<<<<<<<<<<<<
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
   plot. <- plot. &  dev.cur() != 1 
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
#       end
#
############################################################################

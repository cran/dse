# The following should not be needed, but seems necessary for R CMD check
#   to work with the bundle  

 .onLoad  <- function(library, section) {
   require("dse1")
   invisible(TRUE)
   }

############################################################################

#    functions for data analysis      <<<<<<<<<<<<<

############################################################################


phasePlots <- function(data, max.lag=1, diff=FALSE){
# data is a matrix with a variable in each column.(each column is a
# time series), or an object of class TSdata, inwhich case outputData(data) is used.
# trace plots of data and lagged (difference) data (phase space).
# Non-linearities may show up as a non-linear surface, but this is
#   a projection so, for example, a spherical space would not show up.
#  Some sort of cross-section window would show this but require even more plots.
# A good statistical test would be better!

par(mfcol=c(5,5))  #,mar=c(2.1,4.1,4.1,0.1) )
if (is.TSdata(data)) data <- outputData(data)
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

#    functions for model analysis   <<<<<<<<<<<<<

############################################################################


shockDecomposition <- function(model, horizon=30, shock=rep(1,horizon))
{ m <-nseriesInput(model)
  p <-nseriesOutput(model)

  if (is.TSestModel(model)) 
     {data <- TSdata(model)   # just for input
      model <- TSmodel(model)
     }
  else 
    {if( m > 0 ) data <- TSdata(input=matrix(0,horizon,m))
     else        data <- TSdata(output=matrix(0,horizon,p))
    }
   if (!is.TSmodel(model)) stop("shockDecomposition expecting a TSmodel.")
   model$z0 <- NULL    # zero initial conditions
   par(mfrow = c(p, p) , mar = c(2.1, 4.1,3.1, 0.1) )
   for (i in 1:p)
     {outputData(data) <- matrix(0, horizon, p)
      outputData(data)[,i] <- shock   
      z <- l(model,data)
      tfplot(z, reset.screen = FALSE)
     }
invisible()
}

############################################################################

#    functions for forecasting    <<<<<<<<<<<<<

# Class "featherForecasts" has a forecast path from multiple starting points
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

forecast.TSmodel <- function(obj, data,  horizon=36, conditioning.inputs=NULL,
    conditioning.inputs.forecasts=NULL, percent=NULL, ...)
{#  (... further arguments, currently disregarded)
 # obj must be a TSmodel
 # Calculate (multiple) forecasts from the end of data to a horizon determined
 # either from supplied input data or the argument horizon (more details below).
 # In  the case of a model with no inputs the horizon is determined by the
 #   argument horizon.
 # In the case of models with inputs, on which the forecasts
 #  are conditioned, the argument horizon is ignored (except when percent is
 #  specified) and the actual horizon is determined by the inputs in the 
 #  following way:
 # If inputs are not specified by optional arguments (as below) then the default
 #  will be to use inputData(data). This will be the same as the function l() unless
 #  inputData(data) is longer (after NAs are trimmed from each separately) than
 #  outputData(data).
 # Otherwise, if conditioning.inputs is specified it is used for inputData(data).
 #    It must be a time series matrix or a list of time series matrices each
 #    of which is used in turn as inputData(data). The default above is the same as
 #        forecast(model, trimNA(data), conditioning.inputs=trimNA(inputData(data)) )
 # Otherwise, if conditioning.inputs.forecasts is specified it is appended 
 #   to inputData(data). It must be a time series  
 #   matrix or a list of time series matrices each of which is 
 #   appended to inputData(data) and the concatenation used as conditioning.inputs.
 #   Both conditioning.inputs and conditioning.inputs.forecasts should not be
 #   specified.
 # Otherwise, if percent is specified then conditioning.inputs.forecasts are set
 #    to percent/100 times the value of input corresponding to the last period
 #    of outputData(data) and used for horizon periods. percent can be a vector, 
 #    in which case each value is applied in turn. ie c(90,100,110) would would 
 #    give results for conditioning.input.forecasts 10 percent above and below 
 #    the last value of input.

 # The result is an object of class forecast which is a list with 
 #   elements model, horizon, conditioning.inputs, percent, and forecast.
 #   forecast is a list with TSdata objects as elements, one for each element 
 #   in the list conditioning.inputs.
 
 if(!is.TSmodel(obj)) stop("obj must be a TSmodel in forecast.TSmodel")

 output <- trimNA(outputData(data))
 sampleT <- dim(output)[1]

 if (0==nseriesInput(obj))
     {if (0 != (nseriesInput(data)))
              warning("data has input and model does not take inputs.")
      if (!is.null(conditioning.inputs))
          warning("model does not take inputs. conditioning.inputs ignored.")
      if (!is.null(conditioning.inputs.forecasts))
   warning("model does not take inputs. conditioning.inputs.forecasts ignored.")
      if (!is.null(percent))
          warning("model does not take inputs. percent ignored.")
      pr <- l(obj, data, sampleT =sampleT, 
                           predictT=sampleT+horizon)$estimates$pred
      pred <- tfwindow(pr, end=tfend(output), warn=FALSE)
      pr <- tfwindow(pr, start=c(0,1)+tfend(output), warn=FALSE)
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
            {inp <-tframed(rbind(inputData(data),conditioning.inputs.forecasts[[i]]), 
              list(start=tfstart(inputData(data)),
                   frequency=tffrequency(inputData(data))))
             conditioning.inputs <- append(conditioning.inputs, list(inp))
         }  }  
   else if (!is.null(percent))   
        {last.in <- inputData(data)[sampleT,]
         for (i in 1:length(percent) )
           {pol <- t(matrix(last.in*percent[i]/100, length(last.in), horizon))
            inp <-ts(rbind(inputData(data)[seq(sampleT),,drop=FALSE],pol), 
                     start=tfstart(inputData(data)),
                     frequency=tffrequency(inputData(data)))
            conditioning.inputs <- append(conditioning.inputs, list(inp))
        }  }  
   else conditioning.inputs <- trimNA(inputData(data))

   if (is.matrix(conditioning.inputs))
          conditioning.inputs <- list(conditioning.inputs)

   proj <- NULL
   for (policy in  conditioning.inputs)
     {if(!all(tfstart(policy) == tfstart(output)))
            stop("input and output data must have the same starting period (after NAs are removed).")
      pdata <- TSdata(input=policy, output=output)
      predictT <- periods(policy)
      if (sampleT > predictT) 
         stop("input series must be at least as long as output series (after NAs are removed).")
      horizon <- predictT - sampleT
      pr <- l(obj, pdata, sampleT=sampleT, predictT=predictT)$estimates$pred
#    The following lines sometimes cause problems if output is output[...,]
#    See comments in dse2.function.tests
        if (0 == length(proj)) pred <- tfwindow(pr, end=tfend(output), warn=FALSE)
        if(all(tfend(pr)==tfend(output)))
          {pr <- NULL
           warning("Input is not longer than output data. No forecasts produced.") 
          }
        else pr <- tfwindow(pr, start=c(0,1)+tfend(output), warn=FALSE)
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
forecasts <- function(obj)   UseMethod("forecasts")
forecasts.forecast <- function(obj){obj$forecast}


testEqual.forecast <- function(obj1, obj2, fuzz=1e-14)
{# N.B. models are not compared (so equivalent models will compare true)
 # inputs are not compared, so they may be provided differently.
 r <- all(class(obj1) == class(obj2))
 if (r) r <- all(outputData(obj1$data) == outputData(obj2$data))
 if (r) r <- all(obj1$horizon == obj2$horizon)
 if (r) r <- fuzz > max(abs(obj1$pred - obj2$pred))
 if (r) r <- length(obj1$forecast)==length(obj2$forecast)
 for (i in seq(length(obj1$forecast)))
   if (r) r <- fuzz > max(abs(obj1$forecast[[i]] - obj2$forecast[[i]]))
 r
}

tfplot.forecast <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf),
        series = seq(length=nseriesOutput(x$data)),
	Title="Predictions (dotted) and actual data (solid)",
        ylab = seriesNamesOutput(x$data), 
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
{#  (... further arguments, currently disregarded)
 #The default starting point for plots is the start data.
 #The default ending point for plots is the end of forecast.
   if (is.null(x$forecast[[1]]))
      stop("Object to be plotted contains no forecast")
   output <- trimNA(outputData(x$data))
   Ngraphs <- min(length(series), graphs.per.page)
   if(reset.screen) {
      old.par <- par(mfcol = c(Ngraphs, 1), mar= mar, no.readonly=TRUE) #c(5.1,6.1,4.1,2.1))
      on.exit(par(old.par)) }
   N <- length(x$forecast)
   H <- 0
   for (t in 1:N) H <- max(H, dim(x$forecast[[t]])[1])
   tf <-tfExpand(tframe(output), add.end=H)
   for(i in series) 
        {z <- c(output[,i], rep(NA,H))
         for (t in 1:N)
            {zz <- c(rep(NA,periods(output)),x$forecast[[t]][,i],
                     rep(NA,H-dim(x$forecast[[t]])[1]))
             zz[periods(output) ] <- output[periods(output), i] #so line joins last data to first forecast
             z <- cbind(z,zz)
            }
         tframe(z) <- tf
         tfOnePlot(z, ylab = ylab[i],start=start, end=end)
         if(!is.null(Title) && (i == series[1]) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
        }
   invisible()
}

seriesNamesOutput.forecast <- function(x)
   {m <- seriesNamesOutput(x$model)
    d <- seriesNamesOutput(x$data)
    if(!all(m == d))
       warning("data and model names do not correspond. Model names returned.")
    m
   }

seriesNamesInput.forecast <- function(x)
   {m <- seriesNamesInput(x$model)
    d <- seriesNamesInput(x$data)
    if(!all(m == d))
       warning("data and model names do not correspond. Model names returned.")
    m
   }

############################################################################

#    methods for featherForecasts        <<<<<<<<<<<<<

############################################################################

is.featherForecasts <- function(obj) inherits(obj, "featherForecasts")

nseries.featherForecasts <- function(x) nseries(x$featherForecasts[[1]])

seriesNamesOutput.featherForecasts <- function(x) seriesNamesOutput(x$data)
 seriesNamesInput.featherForecasts <- function(x)  seriesNamesInput(x$data)

featherForecasts <- function(obj, ...) UseMethod("featherForecasts")

featherForecasts.TSestModel <- function(obj, data=NULL, ...)
     {if (is.null(data)) data <- TSdata(obj)
      featherForecasts(TSmodel(obj), data, ...)}

featherForecasts.TSdata <- function(obj, model, ...)
     {featherForecasts(model, obj, ...)}

featherForecasts.TSmodel <- function(obj, data, horizon=36,
             from.periods =NULL, ...)
  {data <- freeze(data)
   if (is.null(from.periods))
     {if(0 == nseriesOutput(data)) from.periods <-
             10*seq(floor(periods(data)/10))
      else from.periods <-
             10*seq(floor(min(periods(data), periodsInput(data)-horizon)/10))
     }
   # periods.TSPADIdata returns NA rather than fetching data. Note:Previously freeze
   #  was not done above and pr below just left NULL for TSPADIdata, so some
   #  of this can be cleaned out.
   if ((!is.na(periods(data))) && (max(from.periods) > periods(data) ))
     stop("from.periods cannot exceed available output data.")
   if (0 != (nseriesInput(data)))
     if ((!is.na(periodsInput(data))) && 
        ((max(from.periods)+horizon) > periodsInput(data) ))
       stop("forecasts cannot exceed available input data.")
   shf <- startShift(obj,data,y0=NULL)  # ? y0=y0)
   proj <- NULL
   for (sampleT in from.periods)
     {pr <- l(obj, data, sampleT=sampleT, 
              predictT=sampleT+horizon, result="pred", ...)
      pr[1:(sampleT-1),] <- NA
      # make period before prediction = data so graphics look correct.
      #following if is kludge to prevent retrieving data
      pr[sampleT,] <- outputData(data)[sampleT+shf$shift*shf$lags,]
      proj <- append(proj, list(pr))
     }
   # names are available from  data or obj (model)
   invisible(classed(list(model=obj, # featherForecasts constructor
                data=data, from.periods=from.periods, 
                horizon=horizon, featherForecasts=proj), "featherForecasts"))
}


forecasts.featherForecasts <- function(obj){obj$featherForecasts}

tfplot.featherForecasts <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), 
   series=seq(nseries(x)), 
   Title="Predictions (dotted) and actual data (solid)", 
   ylab=seriesNamesOutput(x),
   graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
{#  (... further arguments, currently disregarded)
 #The default starting point for plots is the start of data.
 #The default ending point for plots is the end of forecasts.
   freq <- tffrequency(x$data)
   names <- seriesNamesOutput(x)
   if(is.null(names)) names <- paste("output", series)
   if (is.null(start)) start <- tfstart(outputData(x$data))
   if (is.null(end))   end   <- addDate(tfend(outputData(x$data)),
                            max(x$horizon), tffrequency(outputData(x$data)))
   if (!is.numeric(series)) series <- match(series, names)
   if(reset.screen) 
     {Ngraphs <- length(series)
      Ngraphs <- min(Ngraphs, graphs.per.page)
      old.par <- par(mfcol = c(Ngraphs, 1), mar= mar, no.readonly=TRUE) #c(5.1,6.1,4.1,2.1)) 
      on.exit(par(old.par))
     }
   # if below is a kludge to skip getting TSPADI data.  
   for(i in series) 
        {#if (is.TSPADIdata(x$data)) 
         #  {zz <- NULL # kludge
         #   ltys <- rep(2,length(x$from.periods))
         #  }
         #else 
           {zz <- tfwindow(outputData(x$data,series=i), start=start,warn=FALSE)
            ltys <- c(1,rep(2,length(x$from.periods)))
           }
         for (t in 1:length(x$from.periods))
            zz <- tbind(zz, selectSeries(x$featherForecasts[[t]], i))
         tfOnePlot(zz,start=start,end=end, ylab=ylab[i], lty=ltys)
         if(!is.null(Title) && (i == series[1]) && (is.null(options()$PlotTitles)
                || options()$PlotTitles))  title(main = Title)
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
#   The main function in this group is MonteCarloSimulations().
#   The main object class in this group is "simulation"

# The second group are for analysing the convergence of estimators. This
#  is extended to functions of estimates (such as model roots). These
#  functions evaluate a single given estimation method with multiple
#  simulated data sets from a given "true" model.
#  The main function in this group is EstEval().
#  The main object classes in this group are
#    "EstEval"
#     c("rootsEstEval","EstEval")
#     c("TSmodelEstEval","EstEval")
#     c("TSestModelEstEval","EstEval")
#     c("coefEstEval","EstEval")           and
#     c("rootsEstEval","EstEval")


# The third group applies multiple estimation techniques to a given data set.
#  This is primarily a utility for other functions 
#  The main (only) function in this group is estimateModels().
#  It returns an object of class c("estimatedModels")


# The fourth group looks at the forecasts and the covariance of forecasts 
#   for multiple horizons. 
#  The simplest case horizonForecasts() which calculates the forecast for 
#    different horizons at all periods in the sample. This is primarily
#    a utility for calculating the forecast error.
#  The next case is is estimatorsHorizonForecastsWRTdata() which 
#   is an extention of horizonForecasts(). It takes specified data 
#   and estimation techniques and calculates forecasts from the estimated 
#   models.
#  The generic function forecastCov() which considers mulitple  
#   models and calculates the cov of forecasts relative to a given data set.
#   It takes a list of models (+ trend, zero)  and calculates the cov 
#   of predictions. It uses forecastCovSingleModel()
#  The next case, forecastCovEstimatorsWRTdata() uses a list of estimation
#   methods to estimate a list of models and calculate the cov of predictions 
#   relative to one given data set.

#   The next case forecastCovWRTtrue() takes a list of models (+ trend,
#   zero)  and calculates the cov of forecasts relative to data sets 
#   simulated with a true.model.
#   The next case, forecastCovEstimatorsWRTtrue simulates data and 
#   uses a list of estimation methods to estimate a list of models, then
#   calculates the cov of predictions relative to other simulated data set.
#  The main object classes in this group are
#     c("estimatorsHorizonForecastsWRTdata") # ? "horizonForecasts")
#     "horizonForecasts"
#     c("multiModelHorizonForecasts","horizonForecasts")
#     "forecastCov"
#     c("forecastCovWRTdata", "forecastCov")
#     c("forecastCovWRTtrue", "forecastCov")
#     c("forecastCovEstimatorsWRTdata",  "forecastCov")
#     c("forecastCovEstimatorsWRTtrue",  "forecastCov")


# The fifth group are some experimental estimation techniques.

############################################################################
#
#       utilities  <<<<<<<<<<
#
############################################################################


minimumStartupLag <- function(model)UseMethod("minimumStartupLag")
minimumStartupLag.TSestModel <- function(model)  
    {minimumStartupLag(model$model)}

minimumStartupLag.ARMA <- function(model)
  {lag <- dim(model$A)[1] 
   if (!is.null(model$C)) lag <- max(lag, dim(model$C)[1])
   lag
  }

minimumStartupLag.SS <- function(model)  { 1+dim(model$F)[2] }

startShift <- function(model,data, y0=NULL)
 {# there is some redundancy between this and  minimumStartupLag which 
  #   should be cleaned up.
  # This function is used to determine the number of lags (and leads) needed
  #   for a model, and whether the data can be padded with zeros or the start
  #   (and end) have to be shifted within the data. Shifting is indicated if the
  #   model has an element $no.zeros (which would be specified if, for example,
  #   the model takes logrithms of data) or the data come from an external data base.
  # As of Nov. 1995 it is used by l.troll, simulate.troll, monte.carlo.troll and
  #  by featherForecasts and forecastCov to determine lags and whether the
  #  starting point is shifted or zeros prepended to the data.
  # if (is.TSPADIdata(data)) shift <- T requires library dsepadi
  if (inherits(data, "TSPADIdata"))  shift <- TRUE
  else if (!is.null(model$no.zeros)) shift <-model$no.zeros
  else shift <- FALSE

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
#       methods for MonteCarloSimulations  <<<<<<<<<<
#
############################################################################

generateSSmodel <- function(m,n,p, stable=FALSE)
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


MonteCarloSimulations <- function(model, simulation.args=NULL, 
           replications=100, rng=NULL, quiet=FALSE, ...) UseMethod("MonteCarloSimulations")

MonteCarloSimulations.default <- function (model, simulation.args=NULL, 
 		replications=100, rng=NULL, quiet=FALSE, ...){
        #  (... further arguments, currently disregarded)
 	if (is.null(rng)) rng <- setRNG()
	else {
		old.rng <- setRNG(rng)
		on.exit(setRNG(old.rng))
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
	# default does not work on array seriesNames(result) <- seriesNames(r)
        if (length(seriesNames(r)) != dim(result)[2])
         stop("length of names (",length(seriesNames(r)),
	      ") does not match number of series(",dim(result)[2],").")
        attr(result,"seriesNames") <- seriesNames(r)
	result <- tframed(result, tfr)
	invisible(classed(list(simulations = result, model = model, 
        rng = rng,  simulation.args = simulation.args, 
        description = "data generated by MonteCarloSimulations.default"), 
        c("MonteCarloSimulations")))
}


MonteCarloSimulations.TSestModel <- function(model, simulation.args=NULL, 
           replications=100, rng=NULL, quiet=FALSE, ...)
  {if (is.null(simulation.args$sd) & is.null(simulation.args$Cov)) 
     simulation.args$Cov <- model$estimates$cov
   if (is.null(simulation.args$input)) simulation.args$input <- inputData(model)
   MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args, 
           replications=replications, rng=rng, quiet=quiet, ...)
  }

MonteCarloSimulations.EstEval <- function(model, simulation.args=NULL,
            replications=100, rng=getRNG(model),  quiet=FALSE, ...)
     MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args,
         replications=replications, rng=rng,  quiet=quiet, ...)
# this looks like a candidate for NextMethod

MonteCarloSimulations.MonteCarloSimulations <- function(model, simulation.args=NULL,
            replications=100, rng=getRNG(model),  quiet=FALSE, ...)
     MonteCarloSimulations(TSmodel(model), simulation.args=simulation.args,
         replications=replications, rng=rng,  quiet=quiet, ...)
# this looks like a candidate for NextMethod


MonteCarloSimulations.TSmodel <- function(model, simulation.args=NULL,
          replications=100, rng=NULL, quiet=FALSE, ...)
#	  Spawn=if (exists(".SPAWN")) .SPAWN else FALSE, ...)
{ 
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 arglist <- append(list(model), simulation.args)
# if (Spawn)
#  {if (!quiet)cat("Spawning processes to calculate ", replications, " replications.\n")
#   assign("sim.forloop.n", replications, where = 1)
#   assign("sim.forloop.result", list(NULL), where = 1)
# #  assign("sim.forloop.model", model, where = 1)
#   assign("sim.forloop.arglist", arglist, where = 1)
#   on.exit(remove(c("sim.forloop.i", "sim.forloop.n", "sim.forloop.result",
#       "sim.forloop.arglist"),where = 1))
#   For(sim.forloop.i = 1:sim.forloop.n, sim.forloop.result[[sim.forloop.i]] <- 
#       do.call("simulate",  sim.forloop.arglist),
#       first=options(warn=-1), sync = TRUE)
#   result <- array(NA, c(dim(sim.forloop.result[[1]]$output),replications))
#   tfr <- tframe(sim.forloop.result[[1]]$output)
#   for (i in 1:replications) result[,,i] <- sim.forloop.result[[i]]$output
#  }
# else {
   #r <- simulate(model, list(...))$output
   r <- do.call("simulate", arglist)$output
   result <- array(NA, c(dim(r),replications))
   tfr <- tframe(r)
   result[,,1] <- r
   if (1 < replications)
     for (i in 2:replications) 
        result[,,i] <- outputData(do.call("simulate", arglist))
#  }
# default does not work on array seriesNames(result) <- seriesNamesOutput(model)
if (length(seriesNamesOutput(model)) != dim(result)[2])
 stop("length of names (",length(seriesNamesOutput(model)),
      ") does not match number of series(",dim(result)[2],").")
attr(result,"seriesNames") <- seriesNamesOutput(model)

result <- tframed(result, tfr)  # my more general multidimensional ts
invisible(classed( # constructor MonteCarloSimulations
         list(simulations=result, model=model, rng=rng, simulation.args=simulation.args,
              description = "data generated by MonteCarloSimulations.TSmodel"),
   c("MonteCarloSimulations") ))
}

is.MonteCarloSimulations <- function(obj) 
   {inherits(obj,"MonteCarloSimulations")}

print.MonteCarloSimulations <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("Simulation with RNG ", x$rng, " from model:\n")
 print(x$model)
 invisible(x)
}

nseriesOutput.MonteCarloSimulations <- function(x)
   {dim(x$simulations)[2]}

nseriesInput.MonteCarloSimulations <- function(x)
  {nseriesInput(x$simulation.args$data)}

tframe.MonteCarloSimulations <- function(x) tframe(x$simulations)

periods.MonteCarloSimulations <- function(x) periods(tframe(x))

seriesNamesOutput.MonteCarloSimulations <- function(x)
   {dimnames(x$simulations)[[2]]}

seriesNamesInput.MonteCarloSimulations <- function(x)
  {seriesNamesInput(x$simulation.args$data)}

testEqual.MonteCarloSimulations <- function(obj1, obj2, fuzz=1e-16)
 {if (length(obj1$result) != length(obj2$result)) r <- FALSE
  else  r <- all(fuzz > abs(obj1$simulations - obj2$simulations))
  r
 }


summary.MonteCarloSimulations <- function(object,
        series=NULL, periods=1:3, ...)
 {#  (... further arguments, currently disregarded)
  stats <- NULL
  if (!is.null(series))
    {if (dim(object$simulations)[3] <20) 
        warning("SD calculation is not very good with so few simulations.")
     names <- seriesNamesOutput(object)
     if(is.null(names)) names <- seriesNamesOutput(object$model)
     if (!is.numeric(series)) series <-match(series, names)
     names <- names[series]
     mn<-apply(object$simulations[periods,series,,drop=FALSE],c(1,2),mean)
     sd<-apply(object$simulations[periods,series,,drop=FALSE],c(1,2),var)^0.5
     stats <- rbind(mn,sd) 
     dimnames(stats)<- list(c(paste("mean period", periods), 
                          paste("S.D. period",periods)), names)
    }
  classed(list(  # constructor summary.MonteCarloSimulations
     description=object$description,
     sampleT=dim(object$simulations)[1],
     p=      dim(object$simulations)[2], 
     simulations=dim(object$simulations)[3], 
     summary.stats=stats,
     rng=getRNG(object)), 
  "summary.MonteCarloSimulations")
}


print.summary.MonteCarloSimulations <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  cat("Object class MonteCarloSimulations\n")
  cat(x$description, "\n")
  cat("periods=",x$sampleT, "variables=", x$p,"simulations=",x$simulations,"\n")
  cat("rng= ", x$rng, "\n")
  if (!is.null(x$summary.stats))   print(x$summary.stats, digits=digits)
  invisible(x)
 }


tfplot.MonteCarloSimulations <- function(x, 
    tf=tframe(x$simulations), start=tfstart(tf), end=tfend(tf),
    series=seq((dim(x$simulations)[2])), 
    select.simulations=seq(dim(x$simulations)[3]),
    graphs.per.page=5, mar=par()$mar, ...)
  {#  (... further arguments, currently disregarded)
   names <- seriesNames(x$simulations)
   tf.p <- tframe(x$simulations) # actual may be different (but not in default)
   Ngraph <- min(length(series), graphs.per.page)
   old.par <-par(mfcol = c(Ngraph, 1), mar= mar, no.readonly=TRUE) #c(5.1,6.1,4.1,2.1))
   on.exit(par(old.par))
   #zz<- matrix(NA, dim(sim)[1], length(x$simulations))
   if (!is.numeric(series)) series <- match(series, names)
   for(i in series) 
        {zz <- (x$simulations)[,i,select.simulations]
         tframe(zz) <- tf.p
	 seriesNames(zz) <- NULL #otherwise name length does not match in tfwindow.default(zz)
         tfOnePlot(zz,start=start,end=end, ylab=names[i]) #tsplot
         if(i == series[1])  title(main = "Monte Carlo Simulations")
        }
   invisible()
}



distribution.MonteCarloSimulations <- function(obj,
     series=seq(dim(obj$simulations)[2]),
     x.sections=TRUE, periods=1:3, graphs.per.page=5, ...)
  {#  (... further arguments, currently disregarded)
  if (dim(obj$simulations)[3] <20) 
     warning("This is not very good with so few simulations.")
   names <- seriesNamesOutput(obj)
   if(is.null(names)) names <- seriesNamesOutput(obj$model)
   if (!is.numeric(series)) series <- match(series, names)
   names <- names[series]
   Ngraph <- min(length(series), graphs.per.page)
   if (x.sections)
       {data <- obj$simulations[periods, series,,drop=FALSE]
        old.par <-par(mfrow =c(Ngraph, length(periods)),
	              mar=c(5.1,6.1,4.1,2.1), no.readonly=TRUE)
       }
   else 
       {old.par <-par(mfrow =c(Ngraph, 1), mar=c(5.1,6.1,4.1,2.1), no.readonly=TRUE)
        mn <- apply(obj$simulations[, series,,drop=FALSE], c(1,2), mean)
        sd <- apply(obj$simulations[, series,,drop=FALSE], c(1,2), var) ^ 0.5
        plt <- array(c(mn, mn+sd, mn-sd, mn+2*sd, mn-2*sd), c(dim(mn),5))
        tf.p <- tframe(obj$simulations)
      }
   on.exit(par(old.par))
   for (i in 1:length(series)) 
    {if (x.sections)
        {for (j in 1:length(periods)) 
          {if (exists("density")) plot(density(data[j,i,]), # -mean[j,i]
                 type="l",xlab=paste("Period",periods[j]),ylab=names[i])
           else if (exists("ksmooth") & !is.R()) plot(ksmooth(data[j,i,], # -mean[j,i],
                              bandwidth=var(data[j,i,])^0.5, kernel="parzen"),
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
#       methods for EstEval.  <<<<<<<<<<
#
############################################################################

#e.bb.ar.100 <- EstEval( mod2, replications=100, 
#               estimation.args=list(estimation="estVARXar", verbose=F))

#e.bb.ls.over <- EstEval( simple.mod, replications=100, 
#   estimation.args=list(estimation="estVARXls", max.lag=6, verbose=F), 
#   criterion="coef")


EstEval <- function( model, replications=100, rng=NULL, quiet=FALSE, 
                       simulation.args=NULL,
                       estimation=NULL, estimation.args=NULL, 
                       criterion ="coef", criterion.args =NULL) 
#		       Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{
 if(is.null(estimation)) stop("estimation method must be specified.")
 if (is.EstEval(model) | is.MonteCarloSimulations(model))
   {rng  <- getRNG(model)
    model<- TSmodel(model)
   }
  
 truth <- do.call(criterion, append(list(model), criterion.args))

 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
# if (!Spawn) {
    if(!quiet) cat("Calculating ", replications, " estimates.\n")
    result <- vector("list",replications)
    for (i in 1:replications)
       {data <- do.call("simulate", append(list(model), simulation.args))
	m   <-  do.call(estimation, append(list(data),  estimation.args))
        result[[i]]<-do.call(criterion, append(list(m), criterion.args))
       }
#   }
# else {
#     if(!quiet)
#	 cat("Spawning processes to calculate ", replications, " estimates.\n")
#     est.forloop <- function(estimation, estimation.args, model, 
#			      simulation.args, criterion, criterion.args)
#	{data <- do.call("simulate", append(list(model), simulation.args))
#	 m   <-  do.call(estimation, append(list(data),  estimation.args))
#	 do.call(criterion, append(list(m), criterion.args))
#	}
#     assign("est.forloop", est.forloop, where = 1)
#     assign("est.forloop.n", replications, where = 1)
#     assign("est.forloop.result", list(NULL), where = 1)
#     assign("est.forloop.estimation", estimation, where = 1)
#     assign("est.forloop.model", model, where = 1)
#     assign("est.forloop.simulation.args", simulation.args, where = 1)
#     assign("est.forloop.criterion", criterion, where = 1)
#     assign("est.forloop.estimation.args", estimation.args, where = 1)
#     assign("est.forloop.criterion.args", criterion.args, where = 1)
#
#     on.exit(remove(c("est.forloop", "est.forloop.i", "est.forloop.n",
#	  "est.forloop.result",  "est.forloop.estimation","est.forloop.model",
#	  "est.forloop.simulation.args", "est.forloop.criterion", 
#	  "est.forloop.estimation.args","est.forloop.criterion.args"),where = 1))
#
#     For(est.forloop.i = 1:est.forloop.n, est.forloop.result[[est.forloop.i ]]<-
#	est.forloop(est.forloop.estimation, est.forloop.estimation.args, est.forloop.model, 
#	est.forloop.simulation.args, est.forloop.criterion, est.forloop.criterion.args),
#	first=options(warn=-1), sync = TRUE)
#     result<-est.forloop.result
#    }
invisible(classed( # constructor EstEval (EstEvals)
      list(result=result,truth=truth,model=model,
           rng=rng, version=version,
           estimation=estimation, estimation.args=estimation.args,
            criterion=criterion,   criterion.args=criterion.args, 
            simulation.args=simulation.args),
     c(paste(criterion,"EstEval",sep=""), "EstEval")))
}


is.EstEval <- function(obj){inherits(obj,"EstEval")}

testEqual.EstEval <- function(obj1, obj2, fuzz = 0)
 {all(as.character(obj1) == as.character(obj2))}

print.EstEval <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("Estimation evaluation with model:\n")
 print(x$model, digits=digits)
 cat("Evaluation criterion: ",x$criterion, "\n")
 invisible(x)
}

summary.EstEval <- function(object, ...)
 {#  (... further arguments, currently disregarded)
  classed(list( # constructor summary.EstEval
     class=class(object),
     estimation=object$estimation,
     estimation.args= if(!is.list((object$estimation.args)[[1]]))
                           object$estimation.args    else    NULL,
     criterion=object$criterion,
     labels=object$labels,
     criterion.args=object$criterion.args,
     replications=length(object$result), 
     true.model=object$model,
     rng=getRNG(object)), 
  "summary.EstEval")
 }


print.summary.EstEval <- function(x, digits=options()$digits, ...)
{ #  (... further arguments, currently disregarded)
  cat("Object of class: ", x$class, "\n")
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

tfplot.EstEval <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf),
        truth= if(is.TSdata(x$truth)) outputData(x$truth) else x$truth,
        series = seq(length=nseries(truth)),
	Title="Estimated (and true) results",
        ylab = seriesNames(truth), remove.mean=FALSE,
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
 {#  (... further arguments, currently disregarded)
  N<-length(x$result)
  old.par <- par(mfcol = c(min(length(series), graphs.per.page), 1),
           no.readonly=TRUE)
  on.exit(par(old.par))
  tf <- tframe(truth)
  truth <- unclass(truth)
  if (remove.mean) 
    {truth <- sweep(truth,2, colMeans(truth))
     ylab <- paste(ylab, "(mean removed)")
    }
  r <- matrix(NA,periods(truth), N)
  for (j in series) {
  	for (i in 1:N)  r[,i] <- x$result[[i]][,j]
	if (remove.mean) r <- sweep(r,2, colMeans(r))
  	tfOnePlot(tframed(tbind(truth[,j], r), tf), ylab= ylab[j],
	              lty=c(1,rep(2, N)), col=c("black",rep("red",N)),
		      start=start,end=end)
        if(!is.null(Title) && (j == series[1]) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
        }
  invisible()
 }



distribution <- function(obj, ...)UseMethod("distribution")

distribution.TSdata <- function(obj, bandwidth=0.2, 
        select.inputs = seq(length= nseriesInput(obj)),
        select.outputs= seq(length=nseriesOutput(obj)), ...)
  {#  (... further arguments, currently disregarded)
   if (0 !=  nseriesInput(obj))
      distribution( inputData(obj), bandwidth=bandwidth, series=select.inputs)
   if (0 != nseriesOutput(obj))
      distribution(outputData(obj), bandwidth=bandwidth, series=select.inputs)
   invisible(obj)
  }

distribution.default <- function(obj, bandwidth=0.2, series=NULL, ...)
  {#  (... further arguments, currently disregarded)
   # obj should be a ts matrix (perhaps this should be a tf method).
   # If series is NULL then all series are ploted.
   # note that this graphic can be fairly misleading:
   #    distribution(runif(1000))  should be uniform
   names <- seriesNames(obj)
   if (!is.matrix(obj) ) obj <- matrix(obj, length(obj), 1)
   if(!is.null(series)) 
     {obj   <-   obj[,series, drop=FALSE]
      names <- names[series]
     }
   par(mfcol=c(ncol(obj),1))
   for ( i in 1:ncol(obj))
      {if      (exists("density")) rd <- density(obj[,i], bw= bandwidth)
       else if (exists("ksmooth") & !is.R()) rd <- ksmooth(obj[,i], bandwidth=bandwidth) 
       else     stop("Neither ksmooth nor density are available.")
       plot(rd, type="l", ylab="density", ylim=c(0, max(rd$y)), xlab=names[i] )
      }
   invisible()
  }



############################################################################
#
#       methods for rootsEstEval  (EstEval)  <<<<<<<<<<
#
############################################################################

summary.rootsEstEval <- function(object, verbose=TRUE, ...)
{#  (... further arguments, currently disregarded)
  nxt <- if (verbose) NextMethod("summary") else NULL
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
  m <- colSums(r)/N
  cov <- r- t(matrix(object$truth, p, N))
  cov <- (t(Conj(cov)) %*% cov)/(N-1)
  ecov <- r- t(matrix(m, p, N))
  ecov <- (t(Conj(ecov)) %*% ecov)/(N-1)
  classed(list( # constructor summary.summary.rootsEstEval
     nxt=nxt,
     conv=conv,
     true.criterion=object$truth,
     mean=m,
     cov=cov,
     ecov=ecov),
  "summary.summary.rootsEstEval")
 }

print.summary.rootsEstEval <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  if (!is.null(x$nxt)) print(x$nxt)
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


tfplot.rootsEstEval <- function(x, ...)UseMethod("plot.rootsEstEval")

plot.rootsEstEval <- function(x, complex.plane=TRUE, cumulate=TRUE, norm=FALSE,
     bounds=TRUE, transform=NULL, invert=FALSE, Sort=TRUE, ...)
{#  (... further arguments, currently disregarded)
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
    {#plot.roots(x$truth, pch="o")
     plot(x$truth, pch="o")
     for (i in 1:N) addPlotRoots(r[i,], pch="*") 
     addPlotRoots(0, pch="+") # addPlotRoots(0+0i, pch="+")
    }
  else
     {if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      if (cumulate) r <- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
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

roots.rootsEstEval <- function(obj, ...)   {obj}

distribution.rootsEstEval <- function(obj, mod=TRUE, invert=FALSE, Sort=FALSE, 
    bandwidth=0.2, select=NULL, ...)
{#  (... further arguments, currently disregarded)
 # if mod is true the modulus is used, otherwise real and imaginary are separated.
 # if invert is true the reciprical is used.
 # if Sort is true then sort is applied (before cumulate). This is of particular interest
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
      if(!is.null(select)) r <- r[,select, drop=FALSE]
      par(mfcol=c(dim(r)[2],1))
      for ( i in 1:dim(r)[2])
         {if      (exists("ksmooth") & !is.R()) rd <- ksmooth(r[,i], bandwidth=bandwidth) 
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
#       methods for coefEstEval (EstEval)  <<<<<<<<<<
#
############################################################################

summary.coefEstEval <- function(object, verbose=TRUE, ...)
  {classed(summary.rootsEstEval(object, verbose=verbose), "summary.coefEstEval")} # constructor
#  (... further arguments, currently disregarded)


print.summary.coefEstEval <- function(x, digits=options()$digits, ...)
  UseMethod("print.summary.rootsEstEval")
#  (... further arguments, currently disregarded)


tfplot.coefEstEval <- function(x, cumulate=TRUE, norm=FALSE, bounds=TRUE,
        invert=FALSE, Sort=FALSE, graphs.per.page = 5, ...){
      #  (... further arguments, currently disregarded)
      N<-length(x$result)
      n <- 0
      for (i in 1:N) n <- max(n, length((x$result)[[i]]))
      r <- matrix(0,N,n)
      plottrue <- TRUE
      for (i in 1:N) {
	ni <- length((x$result)[[i]])
	r[i,1:ni] <- (x$result)[[i]]
	if (ni != n) plottrue <- FALSE
	}
      if (invert) r <- 1/r
      if(norm)    r <- matrix((rowSums(r^2))^.5, N,1)
      r[is.infinite(r)] <- 0
      if (Sort) r <- t(apply(r,1,sort))
      if (n != length(x$truth)) plottrue <- FALSE
      if (plottrue) {
	true.lines <- c(x$truth)
	if (invert) true.lines <- 1/true.lines
	if(norm) true.lines <-sum(true.lines^2)^.5
	if (Sort) true.lines <- sort(true.lines)
	true.lines <-t(matrix(true.lines, length(true.lines),N))
	if (bounds) {
		z  <- r-true.lines
		#Om <- t(z) %*% z/(nrow(z)-1)
		Om <- crossprod(z, z)/(nrow(z)-1)
		Om <- diag(Om)^.5
		Om <- t(matrix(Om, length(Om), N))
		Om <- Om/matrix((1:N)^.5 , N, ncol(Om))
        }
      }
      if (cumulate) r<- apply(r,2,cumsum)/matrix(1:N,N,ncol(r))
      seriesNames(r) <- paste("parm", seq(ncol(r)))
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

roots.coefEstEval <- function(obj, criterion.args=NULL, ...)
{# extract roots criterion 
  model <- obj$model
  truth <-do.call("roots", append(list(model), criterion.args))
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- setArrays(model)
     r <- append(r, 
           list(do.call("roots", append(list(model), criterion.args))))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"roots"
  obj$criterion.args <-criterion.args
  invisible(classed(obj, c("rootsEstEval","EstEval")))# constructor
}

distribution.coefEstEval <- function(obj,  Sort=FALSE, bandwidth=0.2,
	graphs.per.page=5, ...)
{#  (... further arguments, currently disregarded)
 # if Sort is true then sort is applied (before ave). This is of particular interest
 #   with estimation methods like black.box which may not return parameters
 #   of the same length or in the same order.
      N<-length(obj$result)
      n <- length(obj$truth)
      for (i in 1:N) n <- max(n, length((obj$result)[[i]]))
      if (n == length(obj$truth)) plottrue <- TRUE
      else {
         warning("Number of true parameters does not match number estimated.")
	 plottrue <- FALSE
      }
      r <- matrix(0,N,n)
      for (i in 1:N) r[i,1:length((obj$result)[[i]])] <- (obj$result)[[i]]
      true.lines <- c(obj$truth, rep(0,n-length(obj$truth)))
      if (Sort)
        {r <- t(apply(r,1,sort))
         true.lines <- sort(true.lines)
        }
      xlab <-"parameter "
      old.par <- par(mfcol=c(min(graphs.per.page, ncol(r)),1),
                    mar=c(5.1, 4.1, 4.1, 2.1), no.readonly=TRUE)
      on.exit(par(old.par))
      for ( i in 1:ncol(r))
         {if (is.R())     rd <- density(r[,i], bw=bandwidth)
          else            rd <- ksmooth(r[,i], bandwidth=bandwidth) # Splus
          plot(rd, type="l", ylim=c(0, max(rd$y)),
               ylab="density",  xlab=paste(xlab, i) , main="")
          if (plottrue) lines(rep(true.lines[i],2),c(1,0))
         }
      invisible()
}

TSmodel.coefEstEval <- function(obj, ...)
{# rebuild model from coef
  model <- obj$model
  truth <-TSmodel(model)
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- setArrays(model)
     r <- append(r, list(model))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSmodel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSmodelEstEval","EstEval")))
}

TSestModel.coefEstEval <- function(obj)
{# rebuild ... 
  model <- obj$model
  truth <-l(TSmodel(model), data)   # need to regenerate data
  r <- NULL
  for (m in obj$result)
    {coef(model) <- m
     model <- l( setArrays(model), data)
     r <- append(r, list(model))
    }
  ok <- TRUE
  for (m in obj$result)
     ok <- ok & (length(coef(model)) == length(m))  # not perfect but ...
  if (!ok) warning("Parameters do not all correspond to given true model.")
  obj$result<-r
  obj$truth <-truth
  obj$criterion<-"TSestModel"
  obj$criterion.args <-NULL
  invisible(classed(obj, c("TSestModelEstEval","EstEval")))
}

############################################################################
#
#       methods for TSmodelEstEval (EstEval)  <<<<<<<<<<
#
############################################################################

summary.TSmodelEstEval <- function(object, ...)
 {#  (... further arguments, currently disregarded)
  if (is.null((object$result)[[1]]$converged)) conv <- NULL
  else
    {conv <- rep(NA,length(object$result))
     for (i in 1:length(conv)) conv[i] <- (object$result)[[i]]$converged
    }
 # summary(coef(object))
 #summary(roots(object))  these are slow

  classed(list( # constructor summary.TSmodelEstEval
     class=class(object),
     conv=conv,
     default=summary.default(object)),
  "summary.TSmodelEstEval")
 }

print.summary.TSmodelEstEval <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
  cat("Object of class: ",x$class, "\n")
  if (!is.null(x$conv))
        {if(all(x$conv)) cat("All estimates converged.\n")
         else cat(sum(!x$conv)," estimates did not converge!\n")
        }
  print(x$default)
  invisible(x)
}


coef.TSmodelEstEval <- function(object, criterion.args=NULL, ...)
{#  (... further arguments, currently disregarded)
 # extract parameters from models in the list and 
 #   return a list of class coefEstEval EstEval 
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodelEstEval can provide this argument.
  truth <-coef(object$truth)
  r <- NULL
  for (m in object$result) 
     r <- append(r,list(coef(m)))
  if (! is.null((object$result)[[1]]$converged))
    {rc <- rep(NA,length(object$result))
     for (i in 1:length(rc)) rc[i] <- (object$result)[[i]]$converged
     object$result.conv<-rc
    }
  object$result<-r
  object$truth <-truth
  object$criterion<-"coef"
  object$criterion.args <-criterion.args
  invisible(classed(object, c("coefEstEval","EstEval")))
}

roots.TSmodelEstEval <- function(obj, criterion.args=list( randomize=TRUE), ...)
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
  invisible(classed(obj, c("rootsEstEval","EstEval")))
}


tfplot.TSmodelEstEval <- function(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL, ...){
  #  (... further arguments, currently disregarded)
  # extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSmodelEstEval", sep=""), 
               append(list(x), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
  }

############################################################################
#
#       methods for TSestModelEstEval (EstEval)   <<<<<<<<<<
#
############################################################################

summary.TSestModelEstEval <- function(object, ...)
  {classed(summary.TSmodelEstEval(object), "summary.TSestModelEstEval") }  # constructor 
#  (... further arguments, currently disregarded)

print.summary.TSestModelEstEval <- function(x, digits=options()$digits, ...)
   UseMethod("print.summary.TSmodelEstEval")
#  (... further arguments, currently disregarded)


coef.TSestModelEstEval <- function(object, criterion.args=NULL, ...)
{#  (... further arguments, currently disregarded)
 # extract parameters from models in the list and convergence info.
 #   return a list of class coefEstEval EstEval
 # criterion.args is not used. It is provided only so calls from 
 #   summary.TSmodelEstEval can provide this argument.
  truth <-coef(object$truth)
  r <- NULL
  for (m in object$result) r <- append(r,list(coef(m)))
  rc <- rep(NA,length(object$result))
  for (i in 1:length(rc)) rc[i] <- (object$result)[[i]]$converged
  object$result<-r
  object$result.conv<-rc
  object$truth <-truth
  object$criterion<-"coef"
  object$criterion.args <-criterion.args
  invisible(classed(object, c("coefEstEval","EstEval")))
}

roots.TSestModelEstEval <- function(obj, criterion.args=NULL, ...)
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
  invisible(classed(obj, c("rootsEstEval","EstEval")))
}


tfplot.TSestModelEstEval <- function(x, graph.args=NULL,
                       criterion ="coef", criterion.args=NULL, ...){
  #  (... further arguments, currently disregarded)
  # extract criterion and pass to another method with graph.args
  r <- do.call(paste(criterion,".TSestModelEstEval", sep=""), 
               append(list(x), list(criterion.args=criterion.args)))
  do.call("tfplot", append(list(r), graph.args))
  invisible(r)
  }

############################################################################
#
#       function for generating estimatedModels (and methods).   <<<<<<<<<<
#
############################################################################

estimateModels <- function(data, estimation.sample=NULL,
        trend=FALSE, quiet=FALSE, estimation.methods=NULL)
{# Estimate models from data with methods indicated by estimation.methods. 

  if (!is.null(estimation.sample))
    {# is.integer in the next line does not work 
     if (0 != (estimation.sample %%1))
        stop("estimation.sample must be an integer.")
     if (estimation.sample <= 0)
        stop("estimation.sample must be a positive integer.")
     if (nrow(outputData(data)) < estimation.sample)
        stop("estimation.sample cannot be greater than the sample size.")
     outputData(data) <- outputData(data)[1:estimation.sample,, drop=FALSE]
     if (0 != (nseriesInput(data)))
        inputData(data) <- inputData(data)[1:estimation.sample,, drop=FALSE]
    }
   r <-list(estimation.methods=estimation.methods)
  if (trend) r$trend.coef <- lsfit(1:periods(data), outputData(data))$coef
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
  classed(r, "estimatedModels")
}

is.estimatedModels <- function(obj){ inherits(obj,"estimatedModels") }

testEqual.estimatedModels <- function(obj1, obj2, fuzz = 0)
 {all(as.character(obj1) == as.character(obj2))}

print.estimatedModels <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  cat("Estimated models:\n")
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

summary.estimatedModels <- function(object, ...)
 {#  (... further arguments, currently disregarded)
  if (!is.null(object$trend.coef)) cat("trend coef: ", object$trend.coef, "\n")
  if (!is.null(object$multi.model))
    {estimation.names     <- vector(NA, length(object$multi.model))
     estimation.methods   <- vector(NA, length(object$multi.model))
     estimation.converged <- vector(NA, length(object$multi.model))
     for (j in 1:length(object$multi.model))
        {estimation.names[j] <- names(object$estimation.methods)[j]
         estimation.methods[j] <-    object$estimation.methods[[j]]
         estimation.converged[j] <-  object$multi.model[[j]]$converged
        }
    }
  classed(list(  # constructor summary.estimatedModels
     class=class(object),
     trend.coef=object$trend.coef,
     estimation.names=estimation.names,
     estimation.methods=estimation.methods,
     converged=estimation.converged),
  "summary.estimatedModels")
 }

print.summary.estimatedModels <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
  cat("Object of class: ",x$class, "\n")
  cat("Estimated models:\n")
  if (!is.null(x$trend.coef)) cat("trend coef: ", x$trend.coef, "\n")
  if (!is.null(x$multi.model))
    {for (j in 1:length(x$estimation.methods))
        cat(x$estimation.names[j],  x$estimation.methods[j],
	    x$estimation.converged[j],"\n")
    }
  invisible(x)
}



roots.estimatedModels <- function(obj, digits=options()$digits, mod=FALSE, ...)
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

# Class "horizonForecasts"  has 
#  multiple horizons forecasts calculated from every data point.  
#  This is primarily used for calculating forecast errors at different horizons.

############################################################################

#    methods for horizonForecasts        <<<<<<<<<<<<<

############################################################################

is.horizonForecasts <- function(obj) { inherits(obj,"horizonForecasts") }

testEqual.horizonForecasts <- function(obj1, obj2, fuzz=1e-14)
{# N.B. models are not compared (so equivalent models will compare true)
 r <- all(class(obj1) == class(obj2))
 if (r) r <- testEqual(obj1$data, obj2$data)
 if (r) r <- all(obj1$horizons == obj2$horizons)
 if (r) r <- obj1$discard.before == obj2$discard.before
 if (r) r <- fuzz > max(abs(obj1$horizonForecasts - obj2$horizonForecasts))
 if (r) r <- testEqual(obj1$data, obj2$data)
 r
}

horizonForecasts <- function(obj, ...) UseMethod("horizonForecasts")

horizonForecasts.TSestModel <- function(obj, data=NULL, ...) 
{horizonForecasts(TSmodel(obj), if (is.null(data)) TSdata(obj) else data, ...)
}


horizonForecasts.TSdata <- function(obj, model, ...)
{ horizonForecasts.TSmodel(model, obj, ...)}

horizonForecasts.TSmodel <- function(obj, data, horizons=1:4,
	 discard.before=minimumStartupLag(obj), compiled=.DSEflags()$COMPILED, ...)
{#  (... further arguments, currently disregarded)
 # calculate multiple "horizon"-step ahead forecasts 
 # ie. calculate forecasts but return only those indicated by horizons.
 #     Thus, for example, the result of
 #          horizonForecasts(obj, data horizons=c(1,5))    
 #     would be the one-step ahead and five step ahead forecasts.
 # The result is a list of class horizonForecasts with elements model (a 
 #   TSmodel), data, horizons, discard.before, and horizonForecasts.
 # horizonForecasts is an array with three dimension: 
 #   c(length(horizons), dim(obj$data)).
 # Projections are not calculated before discard.before or after
 #   the end of outputData(data).
 # Each horizon is aligned so that horizonForecasts[h,t,] contains the forecast
 #   for the data point outputData(data)[t,] (from horizon[h] periods prior).
 
   if(!checkConsistentDimensions(obj,data)) stop("dimension error\n")
   if (compiled) proj <- horizonForecastsCompiled(obj, data, 
                           horizons=horizons, discard.before=discard.before)
   else
     {TT <-periods(data)
      proj <- array(NA,c(length(horizons),dim(outputData(data))))
      for (t in discard.before:(TT-min(horizons)) )
        {horizons <- horizons[t <= (TT-horizons)]
         z <- l(obj, data, sampleT=t, predictT=TT)$estimates$pred
         for (h in 1: length(horizons) )
             proj[h,t+horizons[h],] <- z[t+horizons[h],]
        }
     }
   dimnames(proj) <- list(NULL, NULL, seriesNamesOutput(data))
   proj <- list(model=obj, data=data, horizons=horizons, 
                discard.before=discard.before, horizonForecasts=proj)
   invisible(classed(proj, "horizonForecasts" ))
}

horizonForecastsCompiled <- function(obj, data, horizons=1:4,
	  discard.before=minimumStartupLag(obj))
  UseMethod("horizonForecastsCompiled")

horizonForecastsCompiled.ARMA <- function( obj, data, horizons=1:4,
	  discard.before=minimumStartupLag(obj))
{if (discard.before < dim(obj$A)[1] )
       warning(paste("Results may be spurious. discard.before should be set higher than the AR order (=", 
                   dim(obj$A)[1]-1, ")."))
 horizons <- sort(horizons)
  p <- nseriesOutput(data)
  TT <- periods(data)
  proj <- array(0,c(length(horizons),TT,p))
  storage.mode(proj) <- "double"
  m <- nseriesInput(obj)
  if (m==0)
     {C <- array(0,c(1,p,1))    # can't pass 0 length array to compiled
      u <- matrix(0,TT,1)
     }
  else
     {C <-    obj$C
      u <- inputData(data)
      if (discard.before < dim(C)[1] )
        warning(paste("Results may be spurious. discard.before should be set higher than the order of C (=", 
                      dim(C)[1]-1, ")."))
     }
  TREND <- obj$TREND
  if (is.null(obj$TREND)) TREND <- matrix(0, TT,p)
  is  <- max(m,p)

  storage.mode(u)  <- "double"
  storage.mode(outputData(data))  <- "double"
  storage.mode(obj$A)  <- "double"
  storage.mode(obj$B)  <- "double"
  storage.mode(C)  <- "double"
  storage.mode(TREND)  <- "double"
   .Fortran("rmaprj",
                  proj=proj,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  ey= array(double(1),dim(outputData(data))), 
                  as.integer( m), 
                  as.integer( p) ,      
                  as.integer( dim(obj$A)[1]),  # 1+order of A  
                  as.integer( dim(obj$B)[1]),  # 1+order of B  
                  as.integer( dim(C)[1]),  # 1+order of C  
                  as.integer(TT),
                  u,
                  outputData(data),	     
                  obj$A,  
                  obj$B,   
                  C,
                  TREND,
                  as.integer(is),  # scratch array dim
                  matrix(double(1),is,is),  # scratch array
                  matrix(double(1),is,is),  # scratch array
                  double(is),         # scratch array
                  integer(is*is),         # scratch array IPIV
                  DUP=.DSEflags()$DUP,
		  PACKAGE="dse1"
		  )$proj
}

horizonForecastsCompiled.SS <- function( obj, data, horizons=1:4,
	 discard.before=minimumStartupLag(obj))
{ horizons <- sort(horizons)
  p <- nseriesOutput(data)
  TT <- periods(data)
  proj <- array(0,c(length(horizons),TT,p))
  storage.mode(proj) <- "double"
     gain <- is.innov.SS(obj)
     n <- dim(obj$F)[2]
     if (discard.before <= n )
         warning(paste("discard.before should probably be set higher than the state dimension (=", n, ")."))
     if (is.null(obj$G))
       {m<-0
        G<-matrix(0,n,1)       # can't call compiled with 0 length arrays
        u <- matrix(0,TT,1)
       }
     else
       {m <- dim(obj$G)[2]
        G <-obj$G
        u <- inputData(data)[1:periods(data),,drop=FALSE]
       } 
     if (gain)     # K or Q,R can be NUll in obj, which messes up compiled
       {K <-    obj$K
        Q <-    matrix(0,1,1)      #not used
        R <-    matrix(0,1,1)      #not used
       }
     else
       {Q <-    obj$Q
        if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square in compiled
        R <-    obj$R
        K <-    matrix(0,n,p)      # this is used
       }
     if(is.null(obj$z0)) z <-rep(0,n)   # initial state
     else  z <-obj$z0
     if(is.null(obj$P0)) P <- diag(1,n) # initial state tracking error 
     else  P <- obj$P0              # this is not used in innov. objs

     storage.mode(u)  <- "double"
     storage.mode(outputData(data))  <- "double"
     storage.mode(obj$F)  <- "double"
     storage.mode(G)  <- "double"
     storage.mode(obj$H)  <- "double"
     storage.mode(K)  <- "double"
     storage.mode(Q)  <- "double"
     storage.mode(R)  <- "double"
     storage.mode(z)  <- "double"
     storage.mode(P)  <- "double"
     IS <- max(n,p)
     .Fortran("kfprj",
                  proj= proj, 
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  ey= matrix(double(1),TT,p), 
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  as.integer(TT),  
                  u, 
                  outputData(data),  
                  obj$F,   
                  G,	
                  obj$H,  
                  K, 
                  Q,	   
                  R,	 
                  as.integer(gain),
                  z,
                  P, 
	          as.integer(IS),           # scratch arrays for KF, IS
	          matrix(double(1),IS,IS),  #A
	          matrix(double(1),IS,IS),  # AA
	          matrix(double(1),IS,IS),  # PP
	          matrix(double(1),n,n),  # QQ
	          matrix(double(1),p,p),  # RR 
	          rep(double(1),IS),  # Z
	          rep(double(1),IS), # ZZ
	          rep(double(1),IS), # WW		   
                  integer(IS*IS),         # scratch array IPIV
 		  DUP=.DSEflags()$DUP,
		  PACKAGE="dse1"
		  )$proj
}


forecasts.horizonForecasts <- function(obj){obj$horizonForecasts}

tfplot.horizonForecasts <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf),
   series=seq(length=nseriesOutput(x$data)), 
   Title="Predictions (dotted) and actual data (solid)", 
   ylab=seriesNamesOutput(x$data), 
   graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...){
   #  (... further arguments, currently disregarded)
   #If series is not NULL then only indicated variables are plotted
   # if start is null it is set to the beginning of the data.
   # if end is null it is set to the end of the data.
   output <-outputData(x$data)
   if (is.null(start)) start <- tfstart(output)
   if (is.null(end)) end <- tfend(output)
   names <- seriesNames(output)
   if (!is.numeric(series)) series <- match(series, names)
   if(reset.screen) 
     {Ngraphs <- length(series)
      Ngraphs <- min(Ngraphs, graphs.per.page)
      old.par <- par(mfcol = c(Ngraphs, 1), mar= mar, no.readonly=TRUE) #c(5.1,6.1,4.1,2.1)) 
      on.exit(par(old.par))
     }
   tf <- tframe(output)
   for(i in series) 
     {#unclass below in because x$horizonForecasts is not tframed and tbind
      #   complains if the frequencies do not match
      zz<- tframed(tbind(unclass(output)[,i],t((x$horizonForecasts)[,,i])), tf)
      tfOnePlot(zz,start=start, end=end, ylab =ylab[i])
      if(!is.null(Title) && (i == series[1]) && (is.null(options()$PlotTitles)
                || options()$PlotTitles))  title(main = Title)
     }
   invisible()
   }

 
############################################################################
#
#       methods for estimatorsHorizonForecastsWRTdata.   <<<<<<<<<<
#
############################################################################

estimatorsHorizonForecastsWRTdata <- function(data, 
                       estimation.sample=.5, horizons=1:12,quiet=FALSE,
                       estimation.methods=NULL)
{ # estimation.sample indicates the part of the data to use for estimation.
  # If estimation.sample is less than or equal 1.0 it is
  # used to indicate the portion of points to use for estimation. 
  # Otherwise it should be an integer and is used to indicate the number
  # of points from the beginning of the sample to use for estimation. 
  
  if (is.null(estimation.methods)) stop("estimation.methods must be specified.")
  if (estimation.sample <= 1.0 )
     estimation.sample <- as.integer(round(estimation.sample*nrow(outputData(data))))
  r <- list(data=data, estimation.sample =estimation.sample, horizons=horizons,  
            estimation.methods=estimation.methods )

  r$multi.model <- estimateModels(data, estimation.sample=estimation.sample, 
              trend=FALSE,quiet=quiet, 
	      estimation.methods=estimation.methods)$multi.model

  r$horizonForecasts <- vector("list", length(estimation.methods))
  for (j in 1:length(estimation.methods))
    r$horizonForecasts[[j]] <- horizonForecasts(l(r$multi.model[[j]],data),
               horizons=horizons,
	       discard.before=minimumStartupLag(r$multi.model[[j]]))
  classed(r, c("estimatorsHorizonForecastsWRTdata")) #? "horizonForecasts")
}


############################################################################
#
#       methods for forecastCov.   (including multiple models)<<<<<<<<<<
#
############################################################################

horizonForecasts.forecastCov <- function(obj,horizons=NULL,
 discard.before=NULL, ...)
{#  (... further arguments, currently disregarded)
 # Calculate forecasts of an object for which cov has been calculated.
 # In a sense this is a step backward, but is sometimes useful to look at
 # forecasts after methods have been analysed on the basis of cov. 
 if(is.null(horizons))       horizons <- obj$horizons
 if(is.null(discard.before)) discard.before <- obj$discard.before
 if (!is.null(obj$model))
   {proj <- horizonForecasts.TSmodel(obj$model, obj$data, horizons=horizons, 
                       discard.before=discard.before)
    class(proj) <- "horizonForecasts"
   }
 else if (!is.null(obj$multi.model))
   {proj <-vector("list", length(obj$multi.model))
    for (i in seq(length(obj$multi.model)))
      proj[[i]] <-horizonForecasts.TSmodel(
             (obj$multi.model)[[i]], obj$data, 
             horizons=horizons, discard.before=discard.before)
    class(proj) <- c("multiModelHorizonForecasts","horizonForecasts")
   }
 else  stop("Object does not include a model.\n")
 invisible(proj)
}

tfplot.multiModelHorizonForecasts <- function(x, 
        tf=NULL, start=tfstart(tf), end=tfend(tf), series=NULL, ...){
  #  (... further arguments, currently disregarded)
  for (i in seq(length(x)))
    {tfplot(x[[i]], start=start, end=end, series=series)
     #cat("press return to continue>");key<-dsescan(what="");cat("\n")
     cat("press Enter to continue>"); key<- readLines(n=1)
    }
  invisible()
  }

TSmodel.forecastCov <- function(obj, select=1, ...)
  {if (is.null(obj$multi.model)) NULL else obj$multi.model[[select]]}

TSdata.forecastCov <- function(data, ...) {data$data}

forecastCov <- function(obj, ..., data=NULL, horizons=1:12, discard.before=NULL,
       zero=FALSE, trend=FALSE, estimation.sample= NULL, compiled=.DSEflags()$COMPILED) UseMethod("forecastCov")
  # Use model and data to construct the cov of predictions at horizons.
   # Discard predictions before (but not including) discard.before to remove 
   #    initial condition problems or do out-of-sample analysis.
   #  obj can be a TSestModel or 
   #        a TSmodel, in which case the second arg must be TSdata, or
   #          TSdata,  in which case the second arg must be a TSmodel.
   
  

forecastCov.TSdata <- function(obj, ..., data=NULL, 
   horizons=1:12, discard.before=1, zero=FALSE, trend=FALSE,
   estimation.sample= NULL, compiled=.DSEflags()$COMPILED)
{#  (... further arguments, currently disregarded). Also included 
 #  as in generic, but currently disregarded: zero, trend, estimation.sample,
 # Use pred$output as the predictions of data and calculate forecastCov
 # This is mainly useful for a fixed prediction like zero or trend.
 # The calculation is dominated by sample effects: more points are
 #  dropped from the end for longer prediction horizons; the trend
 #  predictions are better for the first few periods.
 # With very large samples the result should be nearly constant for 
 # all horizons.
 # The default discard.before=1 should work ok for data, but is not 
 #    consistent with the value for model forecasts. When this routine is
 #    called by other functions the value will usually be overridden.
   pred <- obj
   horizons <- sort(horizons)
   p <- nseriesOutput(data)
   TT  <- periods(data)
   cov <- array(0,c(length(horizons), p,p))
   N <- rep(0,length(horizons))   # the sample size used at each horizon
   err <- pred$output - outputData(data)
   if (compiled)
     {storage.mode(cov) <-"double"
      storage.mode(err) <-"double"
      r <- .Fortran("datepr",
                  forecastCov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  as.integer(p), 
                  predictT=as.integer(TT), 
                  err, 
		  DUP=.DSEflags()$DUP,
		  PACKAGE="dse1"
		  ) [c("forecastCov","sample.size")]
     }
   else
     {for (t in discard.before:(TT-horizons[1]+1))
        {h <- t-1+horizons[(t-1+horizons) <= TT]
         e <- err[h,,drop=FALSE]
         for (k in 1:length(h))
           {N[k] <- N[k]+1
            cov[k,,] <- cov[k,,]*((N[k]-1)/N[k]) + e[k,] %o% e[k,]/N[k] 
           }
        }
       r <- list( forecastCov=cov, sample.size=N)
     }
  dimnames(r$forecastCov) <- list(paste("horizon",as.character(horizons)),NULL,NULL)
  r$forecastCov <- list(r$forecastCov)
  r <- append(r, list(pred=pred, data=data, model=NULL, horizons=horizons, 
                      discard.before=discard.before))
  classed(r, "forecastCov")
}

forecastCov.TSestModel <- function(obj, ..., data=TSdata(obj), 
  horizons=1:12, discard.before=NULL, zero=FALSE, trend=FALSE, 
  estimation.sample= periods(TSdata(obj)), compiled=.DSEflags()$COMPILED)
 {forecastCov(TSmodel(obj), ..., data=data,
     horizons=horizons, discard.before=discard.before, zero=zero, trend=trend,
     estimation.sample=estimation.sample, compiled=compiled)}

forecastCov.TSmodel <- function(obj, ..., data=NULL, 
       horizons=1:12, discard.before=NULL, zero=FALSE, trend=FALSE, 
       estimation.sample= periods(data), compiled=.DSEflags()$COMPILED)
 {if (is.null(data)) stop("data= must be supplied.")
  model.list <- list(obj, ...)
  r <- list(data=data, horizons=horizons, discard.before=discard.before)
  if (is.TSmodel(model.list)) model.list <- list(model.list)
  r$forecastCov <-vector("list", length(model.list))
  i <-0  
  for (model in model.list)
      {i <- i+1
       if (is.null(discard.before))
             rn <-  forecastCovSingleModel(TSmodel(model), data, 
                           horizons=horizons, compiled=compiled)
       else  rn <-  forecastCovSingleModel(TSmodel(model), data, 
                           horizons=horizons, discard.before=discard.before, 
                           compiled=compiled)
       #  $ in the following causes problems for some reason
       r$forecastCov[[i]] <- rn$forecastCov
       r$sample.size   <- rn$sample.size
      }
  if (trend)
     {y <- outputData(data)[1:estimation.sample,]
      pred <- cbind(1,1:periods(data)) %*%
                              (lsfit(1:estimation.sample, y)$coef)
      if (is.null(discard.before))
         r$forecastCov.trend <- forecastCov(TSdata(output=pred), data=data,
             horizons=horizons)$forecastCov[[1]]
      else
         r$forecastCov.trend <- forecastCov(TSdata(output=pred), data=data,
           horizons=horizons,discard.before=discard.before)$forecastCov[[1]]
     }
  if (zero)
     {if (is.null(discard.before))
        r$forecastCov.zero <- forecastCov(TSdata(
             output=array(0,dim(outputData(data)))), data=data, 
             horizons=horizons)$forecastCov[[1]]
      else
        r$forecastCov.zero <- forecastCov(TSdata(
	   output=array(0,dim(outputData(data)))), data=data,
           horizons=horizons,discard.before=discard.before)$forecastCov[[1]]
     }
  r$multi.model <- model.list
  classed(r, c("forecastCovWRTdata", "forecastCov"))
}

is.forecastCovWRTdata <- function(obj){inherits(obj,"forecastCovWRTdata")}

forecastCovSingleModel <- function( model, data=NULL, horizons=1:12, 
          discard.before=minimumStartupLag(model), compiled=.DSEflags()$COMPILED)
{ if(!checkConsistentDimensions(model,data)) stop("dimension error.")
  if (discard.before < 1) stop("discard.before cannot be less than 1.")
  horizons <- sort(horizons)
  names <- seriesNames(data)$output
  if (compiled) 
     r <- forecastCovCompiled(model, data, horizons=horizons,
     		 discard.before=discard.before)
  else
    { p <- nseriesOutput(data)
      shf <- startShift(model,data) #,y0=y0)
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
         e <- pred[h,,drop=FALSE]- outputData(data)[h,,drop=FALSE]
         for (k in 1:length(h))
           {N[k] <- N[k]+1
            cov[k,,] <- cov[k,,]*((N[k]-1)/N[k]) + e[k,] %o% e[k,]/N[k] 
           }
        }
       r <- list( forecastCov=cov, sample.size=N)
     }
  dimnames(r$forecastCov) <- list(paste("horizon",as.character(horizons)),names,names)
#  old:
# The following puts the cov in a sub list. This seems unnecessary for a single
#   cov, but means the same structure can be used with multiple model covs.
#  r$forecastCov <- list(r$forecastCov)
#  r <- append(r, list(model=model, data=data, horizons=horizons, 
#                     discard.before=discard.before))
#  class(r) <- "forecastCov"
 r
}

forecastCovCompiled <- function(model, data, horizons=1:12 ,
    discard.before=minimumStartupLag(model)) 
   if (exists(paste("forecastCovCompiled.", class(model)[1], sep="")))
      UseMethod("forecastCovCompiled") else
      stop("compiled code for this model class is not available. Try forecastCov( ...,compiled=F)")

forecastCovCompiled.ARMA <- function(model, data, horizons=1:12 ,
    discard.before=minimumStartupLag(model))
{ if (discard.before < dim(model$A)[1] )
       warning(paste("Results may be spurious. discard.before should be set higher than the AR order (=",
                    dim(model$A)[1]-1, ")."))
  horizons <- sort(horizons)
  p <- nseriesOutput(data)
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
      u <- inputData(data)
      if (discard.before < dim(C)[1] )
        warning(paste("Results may be spurious. discard.before should be set higher than the order of C (=", dim(C)[1]-1, ")."))
     }
  TREND <- model$TREND
  if (is.null(model$TREND)) TREND <- matrix(0,TT,p)
  storage.mode(cov) <-"double"
  is  <- max(m,p)
  storage.mode(u)  <- "double"
  storage.mode(outputData(data))  <- "double"
  storage.mode(model$A)  <- "double"
  storage.mode(model$B)  <- "double"
  storage.mode(C)  <- "double"
  storage.mode(TREND)  <- "double"
  .Fortran("rmaepr",
                  forecastCov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  pred= array(double(1),dim(outputData(data))), 
                  as.integer( m), 
                  as.integer( p) ,      
                  as.integer( dim(model$A)[1]),  # 1+order of A  
                  as.integer( dim(model$B)[1]),  # 1+order of B  
                  as.integer( dim(C)[1]),  # 1+order of C  
                  predictT=as.integer(TT),
                  as.integer(nrow(outputData(data))), 
                  u,
                  outputData(data),	     
                  model$A,  
                  model$B,   
                  C,
                  TREND,
                  as.integer(is),  # scratch array dim
                  matrix(double(1),is,is),  # scratch array
                  matrix(double(1),is,is),  # scratch array
                  double(is),         # scratch array
                  integer(is*is),         # scratch array IPIV
                  DUP=.DSEflags()$DUP,
		  PACKAGE="dse1"
		  )[c("forecastCov","sample.size")]
}

forecastCovCompiled.innov <- function(model, data, horizons=1:12 ,
    discard.before=minimumStartupLag(model))
     NextMethod("forecastCovCompiled")
#  {forecastCovCompiled.SS()}

forecastCovCompiled.nonInnov <- function(model, data, horizons=1:12 ,
    discard.before=minimumStartupLag(model)) 
     NextMethod("forecastCovCompiled")
#  {forecastCovCompiled.SS()}

forecastCovCompiled.SS <- function(model, data, horizons=1:12 ,
    discard.before=minimumStartupLag(model))
{ horizons <- sort(horizons)
  p <- nseriesOutput(data)
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
        u <- inputData(data)
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
     storage.mode(u) <-"double"
     storage.mode(outputData(data)) <-"double"
     storage.mode(model$F) <-"double"
     storage.mode(G) <-"double"
     storage.mode(model$H) <-"double"
     storage.mode(K) <-"double"
     storage.mode(Q) <-"double"
     storage.mode(R) <-"double"
     storage.mode(z) <-"double"
     storage.mode(P) <-"double"
     IS <- max(n,p)
     .Fortran("kfepr",
                  forecastCov=cov,    
                  as.integer(discard.before), 
                  as.integer(horizons), 
                  as.integer(length(horizons)), 
                  sample.size=as.integer(rep(0, length(horizons))),
                  pred= array(double(1),dim(outputData(data))), 
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  predictT=as.integer(TT), 
                  as.integer(nrow(outputData(data))),  
                  u, 
                  outputData(data),  
                  model$F,   
                  G,	
                  model$H,  
                  K, 
                  Q,	   
                  R,	 
                  as.integer(gain),
                  z,
                  P, 
	          as.integer(IS),           # scratch arrays for KF, IS
	          matrix(double(1),IS,IS),  #A
	          matrix(double(1),IS,IS),  # AA
	          matrix(double(1),IS,IS),  # PP
	          matrix(double(1),n,n),  # QQ
	          matrix(double(1),p,p),  # RR 
	          rep(double(1),IS),  # Z
	          rep(double(1),IS), # ZZ
	          rep(double(1),IS), # WW		   
                  integer(IS*IS),         # scratch array IPIV
 		  DUP=.DSEflags()$DUP,
		  PACKAGE="dse1"
		  ) [c("forecastCov","sample.size")]
}

is.forecastCov <- function(obj)
{ inherits(obj,"forecastCov") }

print.forecastCov <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 for (i in 1:dim((x$forecastCov)[[1]])[3]) 
     {cat("   ",dimnames(x$forecastCov)[[1]][i], "\n")
      z <- NULL
      for (j in 1:length(x$forecastCov) )
         z <- tbind(z,  (x$forecastCov)[[j]][,i,i])
      print(z, digits=digits)
     }
 invisible(x)
}


summary.forecastCov <- function(object, horizons=object$horizons, 
    series=seq(nseriesOutput(object$data)), ...)
{#  (... further arguments, currently disregarded)
 names <- seriesNamesOutput(object$data)
 if(!is.numeric(series)) series <- match(series, names)
 names <- names[series]
 descriptions  <-  vector(NA,    length(object$multi.model))
 summary.stats <-  vector("list",length(object$multi.model))

 for (i in seq(length(summary.stats)))
   {descriptions[i] <- object$multi.model[[i]]$description
    z <- NULL
    for (h in seq(length(horizons))) z <- rbind(z,
             diag(object$forecastCov[[i]][h,series,series])^0.5)
    dimnames(z) <- list(paste("S.D.horizon", horizons), names)
    summary.stats[[i]] <- z
   }
  classed(list(  # constructor summary.forecastCov
     class=class(object),
     horizons=length(object$horizons),
     models=length(object$multi.model),
     seriesNamesOutput=seriesNamesOutput(object$data),
     descriptions=descriptions,
     names=names,
     summary.stats=summary.stats,
     nxt=NextMethod("summary")),
  "summary.forecastCov")
 }



print.summary.forecastCov <- function(x, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  cat("class: ", x$class,"   ")
  cat( length(x$horizons), " horizons\n")
  cat(length(x$models),"models\n")
  for (i in seq(length(x$models)))
    {cat("Model", i, x$description[i],  "\n")
     cat("   variable", x$names,"\n")
     print(x$summary.stats[[i]], digits = digits)   
    }
  cat("\n")
  print(x$nxt)
  invisible(x)
 }





testEqual.forecastCov <- function(obj1, obj2, fuzz=1e-14)
{if (is.null(obj1$rng)) ok <- TRUE
 else ok <- testEqual(obj1$rng , obj2$rng)
 if (ok & !is.null(obj1$forecastCov.true) )
  {if (is.null(obj2$forecastCov.true)) ok <-  FALSE
   ok <- fuzz > max(abs(obj1$forecastCov.true-obj2$forecastCov.true))
  }
 if (ok & !is.null(obj1$forecastCov.zero)) 
  {if (is.null(obj2$forecastCov.zero)) ok <- FALSE
   else ok <- fuzz > max(abs(obj1$forecastCov.zero-obj2$forecastCov.zero))
  }
 if (ok & !is.null(obj1$forecastCov.trend)) 
  {if (is.null(obj2$forecastCov.trend)) ok <- FALSE
   else ok <- fuzz > max(abs(obj1$forecastCov.trend-obj2$forecastCov.trend))
  }
 for (i in 1:length(obj1$forecastCov))
   {if (ok & !is.null((obj1$forecastCov)[[i]])) 
         {if (is.null((obj2$forecastCov)[[i]])) ok <- FALSE 
          else ok <- fuzz > 
               max(abs((obj1$forecastCov)[[i]]-(obj2$forecastCov)[[i]]))
      }
   }
 ok
}


totalForecastCov <- function(obj, select=NULL)
{if (is.null(select)) select <-1:dim((obj$forecastCov)[[1]])[2]
 N <- c( dim((obj$forecastCov)[[1]])[1] ,1,1)
 for (j in 1:length(obj$forecastCov) )
   {z <- apply((obj$forecastCov)[[j]],1,diag)
    # $ causes problems
    obj$forecastCov[[j]] <- array(colSums(z[select,]), N)
   }
 if(!is.null(obj$forecastCov.true))
   {z <- apply(obj$forecastCov.true,1,diag)
    obj$forecastCov.true <- array(colSums(z[select,]), N)
   }
 if(!is.null(obj$forecastCov.zero))
   {z <- apply(obj$forecastCov.zero,1,diag)
    obj$forecastCov.zero <- array(colSums(z[select,]), N)
   }
 if(!is.null(obj$forecastCov.trend)) 
   {z <- apply(obj$forecastCov.trend,1,diag)
    obj$forecastCov.trend <- array(colSums(z[select,]), N)
   }
 invisible(obj)
}



tfplot.forecastCov <- function(x, ..., series = 1:dim(x$forecastCov[[1]])[2], 
    select.cov = 1:length(x$forecastCov), select.true = TRUE, 
    select.zero = TRUE, select.trend = TRUE, y.limit = NULL, line.labels = FALSE, 
    lty = NULL, Legend = NULL, Title = NULL, 
    graphs.per.page = 5, mar=par()$mar, reset.screen=TRUE) {
    #  (... further arguments, currently disregarded)
    p <- dim((x$forecastCov)[[1]])[2]
    Ngraph <- 1 + min(length(series), graphs.per.page)
    if (reset.screen) {
       old.par <- par(mfcol = c(Ngraph, 1), mar = mar, no.readonly=TRUE) #c(5.1, 6.1, 4.1, 2.1))
       on.exit(par(old.par)) }
    par(...)
    if (is.null(lty)) 
        lty <- seq(length(select.cov) + 
            (select.true  & !is.null(x$forecastCov.true)) + 
            (select.zero  & !is.null(x$forecastCov.zero)) + 
            (select.trend & !is.null(x$forecastCov.trend)))
    names <- dimnames((x$forecastCov)[[1]])[[2]]
    if (is.null(names)) 
        names <- paste("variable", 1:p)
    for (i in series) {
        z <- matrix(0, length(x$horizons), length(select.cov))
        for (j in 1:length(select.cov))
            z[, j] <- (x$forecastCov)[[select.cov[j]]][, i, i]
        if (select.trend & !is.null(x$forecastCov.trend)) 
            z <- tbind((x$forecastCov.trend)[, i, i], z)
        if (select.zero & !is.null(x$forecastCov.zero)) 
            z <- tbind((x$forecastCov.zero)[, i, i], z)
        if (select.true & !is.null(x$forecastCov.true)) 
            z <- tbind((x$forecastCov.true)[, i, i], z)
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
        if(!is.null(Title) && (i == 1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles))  title(main = Title)
    }
    if (is.null(Legend)) {
        Legend <- paste("prediction covariance", select.cov[show])
        if (!is.null(x$variable.index)) 
           {nm <- c(seriesNamesOutput(x$all.data), 
                  seriesNamesInput(x$all.data))
            for (i in 1:length(select.cov))
                Legend[i] <- paste(Legend[i], "using", 
                   paste(nm[x$variable.index[select.cov[i],]], collapse=" "))
           }
        if (!is.null(x$selection.index)) 
            Legend <- paste(Legend, "=", x$selection.index[show])
        if (select.trend & !is.null(x$forecastCov.trend)) 
            Legend <- c("trend", Legend)
        if (select.zero & !is.null(x$forecastCov.zero)) 
            Legend <- c("zero", Legend)
        if (select.true & !is.null(x$forecastCov.true)) 
            Legend <- c("true", Legend)
    }
    par(mfg = c(Ngraph, 1, Ngraph, 1)) 
#    if (is.R()) 
#       {box(col = 0) #  temp Rbug workaround
#        legend((par()$usr)[1:2], (par()$usr)[3:4], Legend, lty = lty, 
#            col = 1:6, bty = "y")
#       }    else
    legend((par()$usr)[1:2], (par()$usr)[3:4], Legend, lty = lty, 
               col = lty, bty = "y") # this is a bit of an S/R comp. issue.
        #bty is box
    invisible()
    }


############################################################################
#
#   methods for "forecastCovEstimatorsWRTdata", "forecastCov"  <<<<<<<<<<
#                  (multiple estimators,  given data)
#
############################################################################

outOfSample.forecastCovEstimatorsWRTdata <- function(data,
    zero=FALSE, trend=FALSE,
    estimation.sample=.5, horizons=1:12,quiet=FALSE,
    estimation.methods=NULL, compiled=.DSEflags()$COMPILED)
{ # estimation.sample indicates the portion of the data to use for estimation.
  #If estimation.sample is an integer then it is used to indicate the number
  # of points in the sample to use for estimation. If it is a fracton it is
  # used to indicate the portion of points to use for estimation. The remainder
  # of the sample is used for evaluating forecasts.
  
  if (estimation.sample < 1.0 )
     estimation.sample <- as.integer(round(estimation.sample*nrow(outputData(data))))
  discard.before <- 1+estimation.sample
  forecastCovEstimatorsWRTdata(data, estimation.sample, discard.before,
                       horizons=horizons, zero=zero, trend=trend, quiet=quiet,
                       estimation.methods=estimation.methods, compiled=compiled)
}


forecastCovEstimatorsWRTdata <- function(data, estimation.sample=NULL, 
                       compiled=.DSEflags()$COMPILED, discard.before=10,
                       horizons=1:12, zero=FALSE, trend=FALSE,quiet=FALSE,
                       estimation.methods=NULL)
{# Calculate the forecasts cov of models estimated from data with estimation
 #   methods indicated by estimation.methods  (see estimateModels).
 # estimation.sample is an integer indicating the number of points in the
 #     sample to use for estimation. If it is NULL the whole sample is used.
 # discard.before is an integer indicating 1+the number of points in the
 #     beginning of forecasts to discard for calculating covariances.
 # If zero  is T then forecastCov is also calculated for a forecast of zero.
 # If trend is T then forecastCov is also calculated for a forecast of a linear trend.

  r <- list(data=data, estimation.sample =estimation.sample,
            horizons=horizons, discard.before =discard.before, 
            estimation.methods=estimation.methods)
  models <-estimateModels(data, estimation.sample=estimation.sample, 
                       trend=trend,quiet=quiet,
                       estimation.methods=estimation.methods)
  r$multi.model <- models$multi.model
  if (!is.null(estimation.methods))
    {r$forecastCov <- vector("list", length(estimation.methods))
     for (j in 1:length(estimation.methods))
        {rn <-  forecastCovSingleModel(r$multi.model[[j]], data, 
             compiled=compiled, discard.before=discard.before, horizons=horizons)
         r$forecastCov[[j]] <- rn$forecastCov
         r$sample.size   <- rn$sample.size
        }
    }
  if (zero)
     {r$forecastCov.zero <-forecastCov(TSdata(
             output=array(0,dim(outputData(data)))), data=data, 
                discard.before=discard.before, 
                horizons=horizons)$forecastCov[[1]]
     }
  if (trend)
     {pred <- cbind(1,1:periods(data)) %*% models$trend.coef
      r$forecastCov.trend <- forecastCov(TSdata(output=pred), data=data, 
              discard.before=discard.before, 
              horizons=horizons)$forecastCov[[1]]
     }
  classed(r, c("forecastCovEstimatorsWRTdata", "forecastCov")) # constructor
}

is.forecastCovEstimatorsWRTdata <- function(obj)
  {inherits(obj,"forecastCovEstimatorsWRTdata")}

combine.forecastCovEstimatorsWRTdata <- function(e1,e2)
  {if(! testEqual(e1$data, e2$data)) 
       warning("data is not the same. Second set suppressed.")
   if(! all(e1$estimation.sample == e2$estimation.sample)) 
       warning("estimation.sample's are not the same. Second one suppressed.")
   if(! all(e1$horizon == e2$horizon)) 
       stop("horizon's are not the same.")
   if(e1$discard.before != e2$discard.before) 
       warning("discard.before's are not the same. Second one suppressed.")
   e1$forecastCov <- append(e1$forecastCov, e2$forecastCov)
   e1$estimation.methods <- append(e1$estimation.methods, e2$estimation.methods)
# fix   e1$multi.model <- append(e1$multi.model, e2$multi.model)  
   e1
}


extractforecastCov <- function(e,n) UseMethod("extractforecastCov")
  
extractforecastCov.forecastCovEstimatorsWRTdata <- function(e,n)
  {# select indicated forecastCov
   e$forecastCov <- e$forecastCov[[n]]
   e$estimation.methods <- e$estimation.methods[[n]]
   e$multi.model        <- e$multi.model[[n]]  
   e
}

extractforecastCov.forecastCovEstimatorsFromModel <- function(e,n)
  {# select indicated forecastCov
   e$forecastCov <- e$forecastCov[[n]]
   e$estimation.methods <- e$estimation.methods[[n]]
   e$estimatedModels   <- e$estimatedModels[[n]]  
   e
}


tfplot.forecastCovEstimatorsWRTdata <- function(x, 
    series=1:dim(x$forecastCov[[1]])[2], 
    select.cov=1:length(x$forecastCov),
    select.zero=TRUE, select.trend=TRUE,
    graphs.per.page = 5, mar=par()$mar, reset.screen=TRUE, lty=NULL, ...){
 #  (... further arguments, currently disregarded)
 # ... should be arguments to par(). See tfplot.forecastCov for more details.
 Legend<- paste(names(x$estimation.methods), x$estimation.methods)[select.cov]
 if(select.trend & !is.null(x$forecastCov.trend))
       Legend  <- c("trend",Legend)
 if(select.zero  & !is.null(x$forecastCov.zero))
       Legend  <- c( "zero",Legend)
 # NextMethod maybe?
 tfplot(classed(x,"forecastCov"), series=series, 
        select.cov=select.cov, select.true=FALSE,
        select.zero=select.zero, select.trend=select.trend, Legend=Legend, 
        Title="Prediction variance relative to given data.", 
	graphs.per.page = graphs.per.page, mar=mar, reset.screen=reset.screen, 
	lty=lty)
 invisible()
}


############################################################################
#
#   methods for "forecastCovWRTtrue", "forecastCov"  <<<<<<<<<<
#    given true model, evaluate multiple estimation techniques
#    with multiple simulations for estimation and
#         multiple simulations for forecast
#
############################################################################

forecastCovWRTtrue <- function( models, true.model, 
        pred.replications=1, simulation.args=NULL, quiet=FALSE, 
        rng=NULL, compiled=.DSEflags()$COMPILED,
        horizons=1:12, discard.before=10, trend=NULL, zero=NULL, 
	Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{# models should be a list of models
 # The true model is used to generate more
 # data and for each generated data set the forecasts of the 
 # models are evaluated against the simulated data.
 # if trend is not null it is treated as a model output (forecast) and
 # should be the same dimension as a simulation of the models with 
 # simulation.args. If zero is not null a zero forecast is also evaluated.

 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
      
 if (Spawn & (pred.replications > 1))
   {if(!quiet) 
      cat("Spawning processes to calculate ", pred.replications,
            " forecast replications.\n")
    rep.forloop <- function(models, true.model, simulation.args,
                         horizons, discard.before, zero, trend, compiled=.DSEflags()$COMPILED)
      {data<-do.call("simulate",append(list(true.model), simulation.args))
       r <- NULL
       for (j in 1:length(models))
              {r <- c(r, forecastCovSingleModel(models[[j]],data,
                                  compiled=compiled, horizons=horizons, 
                                  discard.before=discard.before)$forecastCov)
              }
         r.true <- forecastCovSingleModel(true.model,data, compiled=compiled,
                               horizons=horizons,
                               discard.before=discard.before)$forecastCov
         if (is.null(trend)) r.trend <- NULL
         else  r.trend <-forecastCov(TSdata(output=trend), data=data,
	     discard.before=discard.before,horizons=horizons)$forecastCov[[1]]
         if(is.null(zero)) r.zero <- NULL
         else r.zero <- forecastCov(TSdata(
              output=array(0,dim(outputData(data)))), data=data, 
              discard.before=discard.before, horizons=horizons)$forecastCov[[1]]
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
       first=options(warn=-1), sync = TRUE)

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
              {fc <- forecastCovSingleModel(models[[j]],data,
                                  compiled=compiled, horizons=horizons, 
                                  discard.before=discard.before)$forecastCov
               if (i == 1) r[[j]] <- fc
               else        r[[j]] <- r[[j]] + fc
              }
         r.true <- r.true+forecastCovSingleModel(TSmodel(true.model),data, 
                               compiled=compiled, horizons=horizons,
                               discard.before=discard.before)$forecastCov
         if (!is.null(trend))
           r.trend <-r.trend+forecastCov(TSdata(output=trend), data=data,
                                       discard.before=discard.before,
                                       horizons=horizons)$forecastCov[[1]]
         if(!is.null(zero))
            r.zero <- r.zero + forecastCov(TSdata(
                output=array(0,dim(outputData(data)))), data=data, 
                discard.before=discard.before, horizons=horizons)$forecastCov[[1]]
        }
      for (j in 1:length(models))  r[[j]] <- r[[j]]/pred.replications
      r.true  <-  r.true/ pred.replications
      if (!is.null(zero)) r.zero  <-  r.zero/ pred.replications
      if (!is.null(trend)) r.trend <-  r.trend/pred.replications
     }
   invisible(classed(  # constructor (forecastCovWRTtrue)
         list(forecastCov=r, forecastCov.true=r.true, 
           forecastCov.zero=r.zero, forecastCov.trend=r.trend,
           multi.model=models,
           rng=rng, version=version,
           pred.replications=pred.replications,
           horizons=horizons, discard.before=discard.before),
        c("forecastCovWRTtrue", "forecastCov")))
}


forecastCovEstimatorsWRTtrue <- function(true.model, rng=NULL,
                       simulation.args=NULL,
                       est.replications=2, pred.replications=2,
                       discard.before=10, horizons=1:12,quiet=FALSE,
                       estimation.methods=NULL, compiled=.DSEflags()$COMPILED)
#		       Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{# Calculate the forecasts cov of models estimated from simulations of 
 # true.model with estimation methods indicated by estimation.methods (see 
 #       estimateModels). 
 # discard.before is an integer indicating 1+the number of points in the
 # beginning of forecasts to discard for calculating forecast covariances.
 # The returned results has element
 #  $forecastCov.true  $forecastCov.zero $forecastCov.trend containing 
 #    covariances averaged over estimation replications and simulation
 #    replications (forecasts will not change but simulated data will).
 #  $forecastCov a list of the same length as estimation.methods with each
 #    element containing covariances averaged over estimation replications 
 #    and simulation replications.
 #  $estimatedModels a list of length est.replications, with each elements as
 #    returned by estimateModels, thus each element has $multi.model as a
 #    subelement containing models for different estimation techniques.  
 #    So, eg.   $estimatedModels[[2]]$multi.model[[1]]  in the result will
 #    be the model from the first estimation technique in the second replication. 

 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 estimatedModels <- vector("list", est.replications)
 for (i in 1:est.replications)
        {data<-do.call("simulate",append(list(true.model), simulation.args))
         models <-estimateModels(data, trend=TRUE,quiet=quiet,
                       estimation.methods=estimation.methods)
         estimatedModels[[i]] <- models
         rn <- forecastCovWRTtrue( models$multi.model, true.model, 
                    pred.replications=pred.replications, zero=TRUE, quiet=quiet,
                    simulation.args=simulation.args, #Spawn=Spawn,
                    horizons=horizons, discard.before=discard.before,
                    trend=cbind(1,1:periods(data)) %*% models$trend.coef,
                    compiled=compiled)
         if (i==1)
           r<-rn[c("forecastCov","forecastCov.true",
                   "forecastCov.zero","forecastCov.trend")]
         else
          {for (j in 1:length(estimation.methods))
              r$forecastCov[[j]] <-  r$forecastCov[[j]]*(i-1)/i + 
                                     rn$forecastCov[[j]]/i
              r$forecastCov.true  <- r$forecastCov.true *(i-1)/i +  
                                     rn$forecastCov.true/i
              r$forecastCov.zero  <- r$forecastCov.zero *(i-1)/i +  
                                     rn$forecastCov.zero/i
              r$forecastCov.trend <- r$forecastCov.trend*(i-1)/i +  
                                     rn$forecastCov.trend/i
          }
        }

  classed(append(r, # constructor forecastCovEstimatorsWRTtrue
       list(true.model=true.model,estimation.methods=estimation.methods,
         estimatedModels=estimatedModels,
         rng=rng, version=version,
         horizons=horizons, 
         discard.before=discard.before, est.replications=est.replications,
         pred.replications=pred.replications, simulation.args=simulation.args)),
      c("forecastCovEstimatorsWRTtrue", "forecastCov") )
}

is.forecastCovEstimatorsWRTtrue <- function(obj)
 {inherits(obj,"forecastCovEstimatorsWRTtrue")}


print.forecastCovEstimatorsWRTtrue <- function(x, digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
cat("forecastCovEstimatorsWRTtrue\n")
 cat("essential data:", x$essential.data, "\n")
 cat("considering:", seriesNamesOutput(x$all.data), 
                      seriesNamesInput(x$all.data), "\n")
 invisible(x)
}

summary.forecastCovEstimatorsWRTtrue <- function(object, digits=options()$digits, ...)
 {#  (... further arguments, currently disregarded)
  conv <- list()
  Ms <- length(object$estimatedModels)
  for (i in seq(Ms)) 
     conv<- append(conv,object$estimatedModels[[i]]$multi.model[[1]]$converged)
  classed(list( # constructor summary.forecastCovEstimatorsWRTtrue
    class(object), 
    horizons=length(object$horizons), 
    Ms=Ms,
    conv=conv,
    nxt=NextMethod("summary")),
  "summary.forecastCovEstimatorsWRTtrue")
 }


print.summary.forecastCovEstimatorsWRTtrue <- function(x, 
      digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("class: ", x[[1]], "   ")
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

roots.forecastCovEstimatorsWRTtrue <- function(obj, digits=options()$digits,
            mod=FALSE, ...)
 {cat("Estimated models:\n")
  if (!is.null(obj$trend.coef)) cat("trend coef: ", obj$trend.coef, "\n")
  if (!is.null(obj$estimatedModels))
    {r <- vector("list",length(obj$estimatedModels)) 
     for (i in 1:length(obj$estimatedModels))
       {cat("estimation ", i, "\n")
        for (j in 1:length(obj$estimatedModels[[i]]$multi.model))
         {r[[i]] <- vector("list",length(obj$estimatedModels[[i]]$multi.model))
          cat("model ", j, "\n")
          r[[i]][[j]] <- roots(obj$estimatedModels[[i]]$multi.model[[j]])
          if (mod) r[[i]][[j]] <- Mod(r[[i]][[j]])
          print(r[[i]][[j]], digits=digits)
          cat("\n")
         }
       }
    }
  invisible(r)
 }


combine.forecastCovEstimatorsWRTtrue <- function(e1,e2)
  {if(! testEqual(e1$true.model, e2$true.model)) 
       warning("true.models are not the same.")
   if(! testEqual(e1$rng == e2$rng)) 
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
   e1$forecastCov <- append(e1$forecastCov, e2$forecastCov)
   e1$estimation.methods <- append(e1$estimation.methods, e2$estimation.methods)
# fix   e1$multi.model <- append(e1$multi.model, e2$multi.model)  
   e1
}

combine.forecastCov <- function(e1,e2)
  {if(! testEqual(e1$model, e2$model)) 
       warning("models are not the same. Second one suppressed.")
   if(! testEqual(e1$data, e2$data)) 
       warning("data is not the same. Second set suppressed.")
   if(! all(e1$sample.size == e2$sample.size))
       warning("sample.sizes are not the same. Second one suppressed.")
   if(! all(e1$horizons == e2$horizons)) 
       stop("horizon's are not the same.")
   if(e1$discard.before != e2$discard.before) 
       warning("discard.before's are not the same. Second one suppressed.")
   e1$forecastCov <- append(e1$forecastCov, e2$forecastCov)
   e1
}

forecastCovReductionsWRTtrue <- function(true.model, rng=NULL,
                       simulation.args=NULL,
                       est.replications=2, pred.replications=2,
                       discard.before=10, horizons=1:12,quiet=FALSE,
                       estimation.methods=NULL,
                       criteria=NULL, compiled=.DSEflags()$COMPILED)
#		       Spawn=if (exists(".SPAWN")) .SPAWN else FALSE)
{
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 information.criteria <- NULL
 for (i in 1:est.replications)
        {data<-do.call("simulate",append(list(true.model), simulation.args))
         models <-estimateModels(data, trend=TRUE,quiet=quiet,
                       estimation.methods=estimation.methods)
         models$multi.model <-MittnikReducedModels(models$multi.model[[1]]) # use only 1
         crit <- NULL
         for (m in models$multi.model) 
            crit<-rbind(crit,informationTestsCalculations(l(m, data)))
         if(!is.null(criteria))
           {addmodels <- vector("list", length(criteria))
            for (i in 1:length(criteria))
              addmodels[[i]] <- models$multi.model[[order(crit[,criteria[i]])[1]]]
            models$multi.model <- append(addmodels, models$multi.model)
           }
         information.criteria <- append(information.criteria, list(crit))
         rn <- forecastCovWRTtrue( models$multi.model, true.model, 
                    pred.replications=pred.replications, zero=TRUE, quiet=quiet,
                    simulation.args=simulation.args, #Spawn=Spawn,
                    horizons=horizons, discard.before=discard.before,
                    trend=cbind(1,1:periods(data)) %*% models$trend.coef,
                    compiled=compiled)
         if (i==1)
           r<-rn[c("forecastCov","forecastCov.true","forecastCov.zero","forecastCov.trend")]
         else
          {for (j in 1:length(models$multi.model))
             r$forecastCov[[j]] <- r$forecastCov[[j]]*(i-1)/i + rn$forecastCov[[j]]/i
           r$forecastCov.true   <- r$forecastCov.true*(i-1)/i + rn$forecastCov.true/i
           r$forecastCov.zero   <- r$forecastCov.zero*(i-1)/i + rn$forecastCov.zero/i
           r$forecastCov.trend <- r$forecastCov.trend*(i-1)/i + rn$forecastCov.trend/i
          }
        }
  classed(append(r, # constructor forecastCovEstimatorsWRTtrue (forecastCovReductionsWRTtrue)
       list(true.model=true.model,
         estimation.methods=c(criteria,estimation.methods),
         rng=rng, version=version,
         horizons=horizons, 
         discard.before=discard.before, est.replications=est.replications,
         pred.replications=pred.replications, simulation.args=simulation.args, 
         information.criteria=information.criteria)),
       c("forecastCovEstimatorsWRTtrue", "forecastCov"))
}



# Functions in the next group are mainly for evaluating the information <<<<<<<<<<<<<
#    content of data series for predicting other series.           <<<<<<<<<<<<<


############################################################################
#
#  methods for "forecastCovEstimatorsWRTdata.subsets", "forecastCov"
#
############################################################################


is.forecastCovEstimatorsWRTdata.subsets <- function(obj) 
{inherits(obj,"forecastCovEstimatorsWRTdata.subsets")}

print.forecastCovEstimatorsWRTdata.subsets <- function(x, 
                 digits=options()$digits, ...)
{#  (... further arguments, currently disregarded)
 cat("forecastCovEstimatorsWRTdata.subsets\n")
 cat("essential data:", x$essential.data, "\n")
 cat("considering:", seriesNamesOutput(x$all.data), 
                      seriesNamesInput(x$all.data), "\n")
 invisible(x)
}

summary.forecastCovEstimatorsWRTdata.subsets <- function(object, ...)
  {#  (... further arguments, currently disregarded)
   classed(list( class(object),  #summary constructor
        horizons=object$horizons, 
        essential.data=object$essential.data,
        output.names=seriesNamesOutput(object$all.data), 
        input.names =seriesNamesInput(object$all.data),
        nxt=NextMethod("summary")), 
    "summary.forecastCovEstimatorsWRTdata.subsets")
  }

print.summary.forecastCovEstimatorsWRTdata.subsets <- function(x,
                 digits=options()$digits, ...)
  {#  (... further arguments, currently disregarded)
   cat("class: ", x[[1]], "   ")
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

genMineData <- function(umodel, ymodel, uinput=NULL, sampleT=100, 
	unoise=NULL, usd=1,ynoise=NULL, ysd=1, rng=NULL)
{if (is.TSestModel(umodel)) umodel <- TSmodel(umodel)
 if (is.TSestModel(ymodel)) ymodel <- TSmodel(ymodel)
 if(!is.TSmodel(umodel)) stop("genMineData expecting a TSmodel.")
 if(!is.TSmodel(ymodel)) stop("genMineData expecting a TSmodel.")

 if (nseriesInput(ymodel) != nseriesOutput(umodel))
   stop("umodel output dimension must equal ymodel input dimension.")
 
 if(is.null(rng)) rng <- setRNG() # returns setting so don't skip if NULL
 else        {old.rng <- setRNG(rng);  on.exit(setRNG(old.rng))  }
 
 input <- outputData(simulate(umodel, input=uinput, sampleT=sampleT, 
                   noise=unoise, sd=usd))
 r <- TSdata(input  = input,
             output = outputData(simulate(ymodel, input=input,
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
{# make a list of univariate models, one for each series in inputData(data)
 #   for use by build.diagonal.model. 
 n <- nseriesInput(all.data)
 multi.models <- vector("list", n)
 for (i in seq(n))
   {data <-trimNA(TSdata(output= inputData(all.data, series=i)))
    multi.models[[i]] <- TSmodel(estVARXls(data, max.lag=max.lag))
   }
 multi.models
}

build.diagonal.model <- function(multi.models)
{# build one diagonal model from a list of models as returned  by 
 # build.input.models. Uses the AR part only. This can be used by genMineData.
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


plot.mineStepwise <- function(x, ...)
  {#  (... further arguments, currently disregarded)
  cases <- length(x$stepwise$rss)
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


mineStepwise <- function(data, essential.data=1,
      method="efroymson", f.crit=2, intercept=TRUE,
      subtract.means=FALSE,  standardize=FALSE, 
      lags.in=6, lags.out=6, trend=FALSE, plot.=TRUE) 
{  data <- freeze(data)
   m <- ncol(inputData(data))
   p <- ncol(outputData(data))
   if(is.null(m))  m <- 0
   N <- nrow(outputData(data))
   if (standardize)
     {svd.cov <- svd(var(outputData(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if (subtract.means)
    {if(m!=0)
       {input.means<-colMeans(inputData(data))
        inputData(data)<-inputData(data)-t(matrix( input.means, m,N))
       }
     output.means <- colMeans(outputData(data))
     outputData(data)  <- outputData(data) - t(matrix(output.means, p,N))
    }
 # The matrix Past is blocks of data:
 #  [ out-1 | out-2 | ... | out-max.lag | in | in-1 | ... | in-max.lag ]
 # so the result M has a corresponding structure.
 # If there is an input variable (m!=0) then it is shifted to give feedthrough 
 #    in one period. If lags.in=lags.out this has the result that a data point
 #    is lost at the beginning of the input series.
   if(m==0)
     {Past <- matrix(NA,N-lags.out, p*lags.out)
      io.indicator <- c(rep(TRUE, p*lags.out))
      v.indicator  <- c(rep(1:p, lags.out)) 
      lag.indicator  <- c(t(matrix(1:lags.out, lags.out,p))) 
      for (i in 0:(lags.out-1)) 
         Past[,(1+p*i):(p*(1+i))] <-outputData(data)[(lags.out-i):(N-i-1),]
      Present <- outputData(data)[(lags.out+1):N, essential.data]
     }
   else 
     {shift <- max(lags.in+1, lags.out) # start pt. for Present
      Past <- matrix(NA,N-shift+1, p*lags.out+m*(1+lags.in))
      io.indicator <- c(rep(TRUE, p*lags.out), rep(FALSE, m*(1+lags.in)))
      v.indicator  <- c(rep(1:p, lags.out), rep(1:m, (1+lags.in))) 
      lag.indicator<- c(t(matrix(1:lags.out, lags.out,  p)), 
                        t(matrix(0:lags.in, (1+lags.in),m))) 
      for (i in 0:(lags.out-1)) 
        Past[,(1+p*i):(p*(1+i))] <-outputData(data)[(shift-1-i):(N-1-i),]
      for (i in 0:lags.in) 
        Past[,(p*lags.out+1+m*i):(p*lags.out+m*(1+i))] <-
                                   inputData(data) [(shift-i):(N-i),]
      Present <- outputData(data)[shift:N, essential.data]
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
             s.output.indicator=s.output.indicator), "mineStepwise")
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



stripMine <- function(all.data, essential.data=1, 
                       estimation.sample=.5, 
                       discard.before=1, horizons=1:12,quiet=FALSE,
                       estimation.methods=NULL,
                       step.size=NULL)
{# Calculate the predictions cov for essential.data of models estimated 
 # with estimation methods indicated by estimation.methods. 
 # estimation.methods is a list with syntax similar to programs
 #  for comparing estimation methods (eg. estimateModels), BUT ONLY 
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
 #      intermediate.result <- stripMine(data, ...)
 #      intermediate.result <- stripMine(intermediate.result)
 #      intermediate.result <- stripMine(intermediate.result)
 #      result <- stripMine(intermediate.result)
 #  This can be done either interactively or in a batch process, but cannot be
 #  done in a function because the database table is not cleared until the top
 #  level expression is complete.
 #  The class of an intermediate result is stripMine.intermediate.result and
 #  the class of the final result is
 #         c("forecastCovEstimatorsWRTdata.subsets", "forecastCov")
 #  If the final result is used in a call to stripMine then it is just 
 #  returned, so extra calls do not cause errors and are very quick.
 #  This is useful when you are too lazy to calculate the exact number of steps.

  if (class(all.data)[1] == "forecastCovEstimatorsWRTdata.subsets")
       {cat("done.\n")
        return(all.data)
       }
  if (class(all.data)[1] == "stripMine.intermediate.result")
    {r <- all.data$forecastCov
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
     p <- nseriesOutput(all.data)
     m <-  nseriesInput(all.data)
     M <- permute(m + p - length(essential.data) )
     # now combine essential data and permutations of non-essential data.
     if (is.null(M))
       {variable.index<-matrix(essential.data,length(essential.data),1)
        warning("essential.data seems to include all series, which does not make sense in call to stripMine.")
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
   {data<-TSdata(output=outputData(all.data, series=variable.index[i,1:p]), 
                  input= inputData(all.data, series=variable.index[i,(p+1):(p+m)]))
    if(0==length(outputData(data))) 
      stop("The variable selection has mistakenly eliminated all output variables.")
    models <-estimateModels(data, estimation.sample=estimation.sample,
                       trend=TRUE,quiet=quiet,
                       estimation.methods=estimation.methods)
    multi.model <- append(multi.model, list(models$multi.model[[1]]))
    rn <- forecastCov( models$multi.model[[1]], data=data, 
                    horizons=horizons, discard.before=discard.before)
    r<- append(r, list(
           rn["forecastCov"][[1]][[1]][,essential.data,essential.data,drop=FALSE]))
   }
  r<-list(forecastCov=r,all.data=all.data, essential.data=essential.data,
         variable.index=variable.index,
         estimation.methods=estimation.methods,
         multi.model=multi.model,
         horizons=horizons, 
         discard.before=discard.before)
  if (end == nrow(variable.index))
    class(r) <- c("forecastCovEstimatorsWRTdata.subsets", "forecastCov")
  else
    {r<-classed(append(r, list(estimation.sample=estimation.sample, #constructor
             quiet=quiet, step.size=step.size, end=end, m=m,p=p)),
           "stripMine.intermediate.result" )
    }
  r
}

# z <-stripMine(eg1.DSE.data.diff, essential.data=1, 
#      estimation.methods= list(estVARXar=list(max.lag=3))) 


minForecastCov <- function(obj, series=1, verbose=TRUE)
  {#obj is an object as returned by stripMine
   #select the min cov for series only!!! at each horizon and print
   # the returned result is a vector indicating the element of forecastCov which
   # was the min at each horizon. It is suitable as an argument to plot eg:
   #     tfplot(obj, select.cov=minForecastCov(obj))
   # The results of this are similar to the default results of 
   #   selectForecastCov() cov info and information about the horizon
   #   where the model is optimal are given.

   N <- length(obj$forecastCov)
   z <- matrix(0,length(obj$horizons),N)
   for (j in 1:N) z[,j]<-obj$forecastCov[[j]][,series,series]
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


selectForecastCov <- function(obj, series=1, 
    select.cov.best=1,
    select.cov.bound=NULL,
    ranked.on.cov.bound=NULL,
    verbose=TRUE)
  {
   N <- length(obj$forecastCov)
   r <- NULL
   if (!is.null(select.cov.bound))
     if (1 == length(select.cov.bound)) 
       select.cov.bound <- rep(select.cov.bound, length(series))
   for (i in 1:length(series)) 
     {z <- matrix(NA,length(obj$horizons),N)
      for (j in 1:N) 
         z[,j]<-obj$forecastCov[[j]][,series[i],series[i]]
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
       {pred[[j]]           <- obj$forecastCov[[r[j] ]]
#       model[[j]]          <- obj$multi.model[[r[j] ]]
       }
   obj$forecastCov <- pred
#  obj$multi.model <- model
   obj$variable.index <- obj$variable.index[r,, drop=FALSE]
   obj$selection.index <- r
   if (verbose)
     {cat("    model  using subset data series (output | input)\n")
      for (j in 1:length(obj$forecastCov))
         cat( j,"   ", r[j],  "   ", 
             obj$variable.index[j,],"\n")
     }
   invisible(obj)
  }


excludeForecastCov <- function(obj, exclude.series=NULL)
  {# exlude results which depend on the indicated series from a 
   #  (forecastCovEstimatorsWRTdata.subsets forecastCov) object.
   if (!is.null(exclude.series))
     {include<- !apply(0 != obj$variable.index[,exclude.series, drop=FALSE], 1,any)
      obj$forecastCov   <- obj$forecastCov[include]
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

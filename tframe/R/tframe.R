
###########################################################################

#  Internal utilities to limit dependence on class/oldclass

if (is.R())
  { tfclass    <- .Alias(class)
   "tfclass<-" <- .Alias(get("class<-"))
  } else
if (is.S())
  { tfclass     <- class
   "tfclass<-" <- function(x, value){ class(x) <- value ; x }
  }

#if (is.S5())
#  { tfclass     <- oldClass
#   "tfclass<-" <- function(x, value){ oldClass(x) <- value ; x }
#  }


classed <- function(x, cls) {tfclass(x) <- cls; x}

# next should use settf to dispatch on the class method for tf not the 
#  class for x. Setting the class of x first really messes up it dispatches 
#   without having the appropriate attributes of the class.
#tfclassed <- function(x, cls, tf) {  tframe(x) <- tf; tfclass(x) <- cls; x}

###########################################################################

# freeze and its default method really belong with tfpadi (dsepadi) but the
#  generic and default are needed more generally than the database interface,
#  so they are included here.

###########################################################################

freeze <- function(data, ...){
#This function allows for the possiblity of data structures which invoke a
#    call to a database.  It is called by functions which actually use data
# eg:    data <- freeze(data)
# in order to take a snapshot from the database and bring it into a 
# structure which can be used for calculations.
    UseMethod("freeze")
} 
 
freeze.default <- function(data)  
  {if ("character"==mode(data)) freeze(tfPADIdata(data, server="ets")) else data} 

###########################################################################

#  Misc. internal utilities

# Use this with "for (i in seq(length=m) )" as m==0 returns NULL and for does no loops
seqN <- function(N) {if (0==length(N)) NULL else if (N<=0) NULL else seq(N)}

###########################################################################



###########################################################################

# tframe classes and methods       <<<<<<<<<<<<

###########################################################################
        

###############################################

#  generic methods  and defaults <<<<<<<<<<<<

################################################

# start, end, frequency, time and window
# are already generic functions in S with a default 
# method which works for vectors and matrices and data with a tsp attribute.

periods <- function(x)
 {# the length in time of the sequence (number of observations)
    UseMethod("periods")
 }

periods.default <- function(x)
  {if (is.array(x)) return(dim(x)[1])
   else return(length(x))
  }

periods.tsp <- periods.default 

#  tfplot and tfprint below provide generic methods for ploting and printing
#  tf time series objects. Plot methods will probably to some processing
#  and eventually call tfplot.default.

tfplot <- function(obj, ...)  UseMethod("tfplot")


tfplot.default <- function(..., xlab=NULL, ylab=NULL,graphs.per.page=5,
                         start.=NULL, end.=NULL,
			 series=seq(nseries(list(...)[[1]])), mar=par()$mar )
 {obj <- list(...)
  d <- obj[[1]]
  if (!is.tframed(d)) UseMethod("plot")
  else
    {names <- series.names(d)
     Ngraphs <- min(length(series), graphs.per.page)
     old.par <- par(mfcol = c(Ngraphs, 1), mar=mar)  
     on.exit(par(old.par))
     for (i in series)
       {z <- matrix(NA, periods(d), length(obj))
        j <- 0
        for (d in obj)
    	  {if (!is.matrix(d)) d <- tframed(as.matrix(d), tframe(d))
	   if(mode(i)=="character") i <- match(i, names)
	   j <- j + 1
	   z[,j] <- select.series(d, series=i)
	  }
	 tfOnePlot(tframed(z, tframe(d), names=names[i]),
	           xlab=xlab, ylab=ylab, start.=start., end.=end.)
	}
    }
  invisible()
 }
 
tfOnePlot <- function(obj, xlab=NULL, ylab=NULL, start.=NULL, end.=NULL, ...)
 {if (!is.tframed(obj)) UseMethod("plot")
  else
    {if (!is.null(start.)) obj <- tfwindow(obj, start. = start.)
     if (!is.null(end.))   obj <- tfwindow(obj, end.   = end.)
     tline <- time(obj)
     if(is.null(xlab)) xlab <- ""
     if(is.null(ylab)) ylab <- paste(series.names(obj), collapse="  ")
     matplot(tline, obj, type="l", xlab=xlab, ylab=ylab, ...)
    }
  invisible()
 }

# Note tfprint prints the data. tframePrint  prints the tframe info. 

tfprint <- function(x, ...)  UseMethod("tfprint")

tfprint.default <- function(x,...)
 {xx <- x
  if(1 == nseries(xx)) xx <- matrix(xx, length(xx), 1)
  dimnames(xx) <- list(format(time(tframe(x))), series.names(x))
  tframe(xx) <- NULL
  series.names(xx) <- NULL
  print(xx, ...)
  invisible(x)
 }



tfwindow <- function(x, ...)  UseMethod("tfwindow")

tfwindow.default <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=T)
  {# this provides a convenient way to support warn and correct for bugs
   # in some versions of window().
   # if (is.null(start.)) start. <- start(x)
   # if (is.null(end.))   end.   <- end(x)
   # With the default warn=T warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   if (is.null(start.) && !is.null(tf)) start. <- start(tf)
   if (is.null(end.)   && !is.null(tf)) end.   <- end(tf)
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   y <- window(x, start=start., end=end.)
   if (is.matrix(x) && !is.matrix(y) )
      y <- tframed(matrix(y, length(y), ncol(x)), tframe(y))
   series.names(y) <- series.names(x)
   y
  }

###############################################

#  tframe  methods   <<<<<<<<<<<<

################################################
is.tframe <- function(tf) inherits(tf, "tframe")
is.tframed <- function(x) inherits(tframe(x), "tframe")

tframe <- function(x) {UseMethod("tframe") } #extract the tframe

tframe.default <- function(x){
	if(is.null(x)) NULL 
	else if (!is.null(attr(x, "tframe"))) attr(x, "tframe") # constructor
	else if (!is.null(tsp(x)))    classed(tsp(x), "tframe") # constructor
	else if(is.array(x) && !is.matrix(x)) 
		classed(tsp(as.ts(seq(dim(x)[1]))), "tframe") # constructor
	else 	classed(tsp(as.ts(x)), "tframe") # constructor
}

# Following switches dispatch to class of tf value rather than class of x
#    except in the case of null assignment.

"tframe<-" <- function(x, value)
 {if(is.null(value)) UseMethod("tframe<-")
  else settf(value,x)
 }

"tframe<-.default" <- function(x, value) {tsp(x) <- value; x}

settf <- function(value, x) {UseMethod("settf") }

settf.default <- function(value, x)
{if (!check.tframeConsistent(value, x))
    stop("time frame value in tframe assignment is not consistent with data.")
 # using tsp is a bit dangerous in R as it sets class "ts" which can clobber
 # other classes. The "pure" idea of tframe is that time frame inheritance
 # follows from the classes of the tframe attribute and classes of the
 # object itself are separate. Instead of setting tsp(x) one would
 #    attr(x, "tframe") <- value
 # This would require start, end, etc to check is.tframed and then dispatch, and
 #  that has not been implemented.
 tsp(x) <- value   
 x
}


tframed <- function(x, ...) {UseMethod("tframed") }

tframed.default <- function(x, tf=NULL, names = NULL) 
{# return x as a tframed object with tframe tf
 # If ts is not a tframe but a list then ts() is attempted. This is not
 #     really the way tframed is suppose to be used, but makes backward 
 #     compatability easier.
    if (!is.null(names))  series.names(x) <-  names
    if (is.null(tf)) tf <- tframe(x)
    if (is.tframe(tf)) x <- settf(tf, x)
    else if (is.list(tf))
       {if( is.null(tf$start) & is.null(tf$end) )
           stop("tf must be a tframe or a list of arguments for ts().")
        x <- do.call("ts", append(list(x), tf))
       }
    x
}



###############################################

#  Generic .tframe methods (these act on the tframe not on the data)

###############################################


tfprint.tframe <- function(x, ...) UseMethod("tframePrint")
tframePrint <- function(x, ...) UseMethod("tframePrint")

tframePrint.default <- function(x, digits=NULL, quote=T, prefix="", ...) 
  {if (! is.tframe(x)) x <- tframe(x)
   invisible(print(unclass(x), quote=quote, prefix=prefix, ...)) }


start.tframe <- function(x)UseMethod("tframeStart")
tframeStart <- function(x)UseMethod("tframeStart")

tframeStart.default <- function(x)
  {if (! is.tframe(x)) x <- tframe(x)
   c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))}


end.tframe <- function(x)UseMethod("tframeEnd")
tframeEnd <- function(x)UseMethod("tframeEnd")

tframeEnd.default <- function(x)
  {if (! is.tframe(x)) x <- tframe(x)
   c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))}


# periods should give the number of data points in the time direction.
periods.tframe <- function(x)UseMethod("tframePeriods")
tframePeriods <- function(x)UseMethod("tframePeriods")

tframePeriods.default <- function(x)
  {if (! is.tframe(x)) x <- tframe(x)
   1+round((x[2]-x[1])*x[3])}


# frequency is less essential and may not always make sense.
frequency.tframe <- function(x)UseMethod("tframeFrequency")
tframeFrequency <- function(x)UseMethod("tframeFrequency")

tframeFrequency.default <- function(x)
  {if (! is.tframe(x)) x <- tframe(x)
   x[3]}


time.tframe <- function(x)UseMethod("tframeTime")
tframeTime <- function(x)UseMethod("tframeTime")

tframeTime.default <- function(x)
  {if (! is.tframe(x)) x <- tframe(x)
   x[1] + (seq(periods(x))-1)/x[3]}


tftruncate.tframe <- function(x, start=NULL, end=NULL)
   UseMethod("tframeTruncate")


tframeTruncate <- function(x, start=NULL, end=NULL) UseMethod("tframeTruncate")
   #NULL means no truncation. 

tfexpand.tframe <- function(x, add.start=0, add.end=0) UseMethod("tframeExpand")

tframeExpand <- function(x, add.start=0, add.end=0)
 UseMethod("tframeExpand")


check.tframeConsistent <- function(tf, x) UseMethod("check.tframeConsistent")

check.tframeConsistent.default <- function(tf, x)
   {tframePeriods(tf) == periods(x)}

test.Equaltframes <- function(tf1, tf2) UseMethod("test.Equaltframes")

test.Equaltframes.default <- function(tf1, tf2) { all(tf1==tf2)}

# Following could be used to do date comparisons like start() < end()

tframeEarliestStartIndex <- function(x, ...)
    UseMethod("tframeEarliestStartIndex")

earliest.tframeStart <- function(x, ...)
    append(list(x),list(...))[[tframeEarliestStartIndex(x, ...)]]

tframeEarliestEndIndex <- function(x, ...)
    UseMethod("tframeEarliestEndIndex")

earliest.tframeEnd <- function(x, ...)
    append(list(x),list(...))[[tframeEarliestEndIndex(x, ...)]]

tframeLatestStartIndex <- function(x, ...)
    UseMethod("tframeLatestStartIndex")

latest.tframeStart <- function(x, ...)
    append(list(x),list(...))[[tframeLatestStartIndex(x, ...)]]

tframeLatestEndIndex <- function(x, ...)
    UseMethod("tframeLatestEndIndex")

latest.tframeEnd <- function(x, ...)
    append(list(x),list(...))[[tframeLatestEndIndex(x, ...)]]


###############################################

#  default .tframe methods   <<<<<<<<<<<<

################################################


tframeTruncate.default <- function(x, start=NULL, end=NULL) 
    {# like window but uses indexes rather than dates
     if (!is.null(end))   x[2] <- x[1] + (end-1)/x[3]
     if (!is.null(start)) x[1] <- x[1] + (start-1)/x[3]
     x
    }

tframeExpand.default <- function(x, add.start=0, add.end=0) 
    {x[2] <- x[2] + add.end/x[3]
     x[1] <- x[1] - add.start/x[3]
     x
    }


tframeEarliestStartIndex.default <- function(x, ...) 
    {r <- 1
     fr <- frequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[1] < args[[r]][1]) r <- i
         }           
     r
    }

tframeEarliestEndIndex.default <- function(x, ...) 
    {r <- 1
     fr <- frequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[2] < args[[r]][2]) r <- i
         }           
     r
    }

tframeLatestStartIndex.default <- function(x, ...) 
    {r <- 1
     fr <- frequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[1] > args[[r]][1]) r <- i
         }           
     r
    }

tframeLatestEndIndex.default <- function(x, ...) 
    {r <- 1
     fr <- frequency(x)
     args <- list(x, ...)
     for (i in seq(length(args)))
         {tf <- args[[i]]
          if (tf[3] != fr) stop("frequencies must be that same.")
          if (tf[2] > args[[r]][2]) r <- i
         }           
     r
    }


###############################################

#   .ts and .tframe.tstframe methods

###############################################

#"tframe.tframe<-.ts" <- function(x, value) {tsp(x) <- value; x}

tframe.ts <- function(x){classed(tsp(x), c("tstframe", "tframe"))} # constructor

"tframe<-.ts" <- function(x, value) {tsp(x) <- value; x}

#settf.default works for .ts

tframeStart.tstframe <- function(x)
   {c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))}

tframeEnd.tstframe <- function(x)
   {c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))}

tframePeriods.tstframe <- function(x)  {1+round((x[2]-x[1])*x[3])}

tframeFrequency.tstframe <- function(x) {x[3]}

tframeTime.tstframe <- function(x) {x[1] + (seq(periods(x))-1)/x[3]}



###############################################

#   .tf and .tframe.tftframe methods

# these provide .ts or tsp like methods for more general arrays with time as the
#   first dimension. It may be possible extend to more general structures.

###############################################

tframe.tf <- function(x) {attr(x, "tf") }

"tframe<-.tf" <- function(x, value)
  {if (!is.null(value)) return(settf.tftframe(value, x))
   else 
     {attr(x, "tf") <- NULL
      if (!is.null(class(x))) class(x) <- class(x)[! class(x) == "tf"]
     }
   x
  }

#tf <- function(x, tfr=NULL, ...){stop("defunct. Use tframed.")}
#tf <- function(x, tfr=NULL, ...)
# {#  tfr should be a tframe or other args are passed to ts.
#  # the class is set to "tftframe" so objects are interchangable R to S
#  if (is.null(tfr)) 
#      {tfr <- tframe(do.call("ts", append(list(seq(periods(x))), list(...))))
#       cls <- class(tfr)
#       cls[cls == "tstframe"] <- "tftframe"
#       class(tfr) <- cls
#      }
#  if (!is.tframe(tfr)) stop("tfr must be a tframe")
#  settf.tftframe(tfr, x)
# }


settf.tftframe <- function(value, x)
{#class(value) <- c("tftframe", "tframe")
 if (!check.tframeConsistent(value, x))
    stop("time frame value in tframe assignment is not consistent with data.")
 attr(x, "tf") <- value
 classed(x, "tf")  # constructor (settf.tftframe)
}


start.tf <- function(x) {start(tframe(x))}
end.tf <- function(x) {end(tframe(x))}
periods.tf <- function(x) {periods(tframe(x))}
frequency.tf <- function(x) {frequency(tframe(x))}
time.tf <- function(x) {time(tframe(x))}

tframeStart.tftframe <- function(x)
   {c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))}

tframeEnd.tftframe <- function(x)
   {c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))}

tframePeriods.tftframe <- function(x)  {1+round((x[2]-x[1])*x[3])}

tframeFrequency.tftframe <- function(x) {x[3]}

tframeTime.tftframe <- function(x) {x[1] + (seq(periods(x))-1)/x[3]}

tfwindow.tf <- function(x, start.=NULL, end.=NULL, warn=T, eps=.Options$ts.eps) 
  {# this needs work
   tfwindow(ts(x, start=start(x), frequency=frequency(x)),
             start.=start., end.=end., warn=warn)
  }

###############################################

#   .TSdata methods (for DSE)

###############################################


"tframe<-.TSdata" <- function(x, value)
 {if (0 != input.dimension(x)) tframe( input.data(x)) <- value
  if (0 !=output.dimension(x)) tframe(output.data(x)) <- value
  x
 }



###############################################

#   .rts and .tframe.rtstframe methods

###############################################



"tframe<-.rts" <- function(x, value) {rts(x) <- value; x}
"tframe.tframe<-.rts" <- function(value, x)
{if (!check.tframeConsistent(value, x))
    stop("time frame value in tframe assignment is not consistent with data.")
 rts(x) <- value
 x
}

###############################################

# add cts its

###############################################


###############################################

#  stamped specific methods   <<<<<<<<<<<<
#  stamped class TS have a date/time stamp associated with each time point
################################################
check.tframeConsistent.stamped <- function(tf, x)
  {periods(x) == periods(tf)}

test.Equaltframes.stamped <- function(tf1, tf2)
  {all(tf1$stamp == tf2$stamp)}

tframePeriods.stamped <- function(x)length(tframe(x))



###############################################



test.equal <- function(obj1, obj2, fuzz = 0) UseMethod("test.equal")

 
test.equal.default <- function(obj1, obj2, fuzz=1e-16) 
  {if      (is.null(obj1)) is.null(obj2)
   else if (is.array(obj1)) test.equal.array(obj1, obj2, fuzz=fuzz)
   else if (is.numeric(obj1)) test.equal.numeric(obj1, obj2, fuzz=fuzz)
   else if (is.list(obj1)) test.equal.list(obj1, obj2, fuzz=fuzz)
   else is.logical(all.equal(obj1, obj2, tolerance=fuzz))
  }

test.equal.array <- function(obj1, obj2, fuzz=1e-16) 
  {if(!is.array(obj2))                     r <-FALSE
   else if (any(dim(obj1) != dim(obj2)))   r <- FALSE
   else if ("character" == mode(obj1))     r <- all(obj1 == obj2)
   else if ("numeric" == mode(obj1))
              r <- test.equal.numeric(obj1, obj2, fuzz=fuzz)
   else stop(paste("matrix of mode ", mode(obj1), " not testable."))
   if (is.na(r))  r <- FALSE
    r
  }

test.equal.matrix <- test.equal.array

test.equal.numeric <- function(obj1, obj2, fuzz=1e-16) 
  {r <- all(is.infinite(obj1) == is.infinite(obj2))
   if (r) 
          {nna <- !is.na(c(obj1))
           r <- fuzz >= max(abs(c(obj1)[nna] - c(obj2)[nna]))
          }
   if (is.na(r))  r <- FALSE
   r
  }

test.equal.list <- function(obj1, obj2, fuzz=1e-16) 
  {r <- length(obj1) == length(obj2)
   if (r) for (i in seq(length(obj1)))
        {if(r) r <- test.equal(obj1[[i]], obj2[[i]], fuzz=fuzz) }
   r
  }

if (!exists("lag")) lag <- function(x, ...) { UseMethod("lag") }

if (!exists("lag.default"))  lag.default <- function(x, ...) {stop("no lag function") }




splice <- function(mat1, mat2, ...) UseMethod("splice")

splice.default <- function(mat1, mat2)
{# splice together 2 time series matrices. If data  is provided in both for
 #  a given period then mat1 takes priority.
 # The result starts at the earlier of mat1 and mat2 and ends at the later.
 # dimnames are taken from mat1.
 # The frequencies should be the same.
 if (is.null(mat1)) return(mat2)
 if (is.null(mat2)) return(mat1)
 freq <- frequency(mat1)
 if (freq != frequency(mat2)) stop("frequencies must be the same.")
 p <- dim(mat1)[2]
 if (p != dim(mat2)[2])   stop("number of series must be the same.")
 fr <- c(freq,1)
 st <- min(fr %*% start(mat1), fr %*% start(mat2))
 strt <- c(st %/% freq, st %% freq)
 en <- max(fr %*% end(mat1), fr%*% end(mat2))
 r1 <-r2 <-tframed(matrix(NA, 1+en-st, p), list(start=strt, frequency=freq))
 r1[c((fr %*% start(mat1))-st) + 1:dim(mat1)[1],] <- mat1
 r2[c((fr %*% start(mat2))-st) + 1:dim(mat2)[1],] <- mat2
 na <- is.na(r1)
 r1[na] <- r2[na] # put mat2 only in na locations of mat1
 dimnames(r1)<-list(round(time(r1),digits=3),dimnames(mat1)[[2]])
 r1 <- tframed(r1, list(start=earliest.start(mat1,mat2), 
                        end =latest.end(mat1,mat2), frequency=freq))
 r1
}


if( !exists("tsmatrix.default"))  
  {if(exists("tsmatrix")) tsmatrix.default <- tsmatrix 
   else tsmatrix.default <- function(x, ...) 
            {tbind(x, ..., pad.start=F, pad.end=F) }
  }

tsmatrix <- function(x, ...)
 {# the default tsmatrix messes up because it gets some time info. (from
  #  start or end) but not tsp info.
  if (is.tframed(x)) tbind(x, ..., pad.start=F, pad.end=F)
  else 
    {#warning("Using tsmatrix which should be defunct. Consider using tbind and tframe methods.")       
     tsmatrix.default(x,  ...)
    }
 }


tftruncate <- function(x, start=NULL, end=NULL)
 {# similar to window but start and end specify periods relative to the 
  #   beginning (eg x[start:end] for a vector).
  #   NULL means no truncation.
  UseMethod("tftruncate")
 }


tftruncate.default <- function(x, start=NULL, end=NULL)
    {tf <- tframeTruncate(tframe(x), start, end)
     if (is.null(start)) start <- 1
     if (is.matrix(x)) 
        {if (is.null(end)) end <- dim(x)[1]
         z <- x[start:end,,drop=F]
        }
     else 
        {if (is.null(end)) end <- length(x)
         z <- x[start:end]
        }
     tframe(z) <- tf
     z
    }

tfexpand <- function(x, add.start=0, add.end=0)
 {# expand (a tframe) by add.start periods on the beginning
  # and add.end periods on the end
  UseMethod("tfexpand")
 }

tfexpand.default <- function(x, add.start=0, add.end=0)
    {tf <- tframeExpand(tframe(x), add.start=add.start, add.end=add.end)
     select.series(tbind(x,time(tf)), series=1)
    }


earliest.start <- function(x, ...)
    start(append(list(x),list(...))[[earliest.startIndex(x, ...)]])

earliest.startIndex <- function(x, ...)
  {if (is.tframe(x)) UseMethod("tframeEarliestStartIndex")
   else 
     {tf <- list(tframe(x))
      for (i in list(...)) tf <- append(tf, list(tframe(i)))
      r <- do.call("tframeEarliestStartIndex", tf)
     }
   r
  }

earliest.end <- function(x, ...)
    end(append(list(x),list(...))[[earliest.endIndex(x, ...)]])

earliest.endIndex <- function(x, ...)
  {if (is.tframe(x)) UseMethod("tframeEarliestEndIndex")
   else 
     {tf <- list(tframe(x))
      for (i in list(...)) tf <- append(tf, list(tframe(i)))
      r <- do.call("tframeEarliestEndIndex", tf)
     }
   r
  }

latest.start <- function(x, ...)
    start(append(list(x),list(...))[[latest.startIndex(x, ...)]])

latest.startIndex <- function(x, ...)
  {if (is.tframe(x)) UseMethod("tframeLatestStartIndex")
   else 
     {tf <- list(tframe(x))
      for (i in list(...)) tf <- append(tf, list(tframe(i)))
      r <- do.call("tframeLatestStartIndex", tf)
     }
   r
  }

latest.end <- function(x, ...)
    end(append(list(x),list(...))[[latest.endIndex(x, ...)]])

latest.endIndex <- function(x, ...)
  {if (is.tframe(x)) UseMethod("tframeLatestEndIndex")
   else 
     {tf <- list(tframe(x))
      for (i in list(...)) tf <- append(tf, list(tframe(i)))
      r <- do.call("tframeLatestEndIndex", tf)
     }
   r
  }



trim.na <- function(x, start. = TRUE, end. = TRUE) UseMethod("trim.na") 

trim.na.default <- function(x, start. = TRUE, end. = TRUE)
{# trim NAs from the ends of a ts matrix or vector.
 # (Observations for all series are dropped in a given period if any 
 #  one contains an NA in that period.)
 # if start.=F then beginning NAs are not trimmed.
 # If end.=F   then ending NAs are not trimmed.
 sample <- ! if (is.matrix(x)) apply(is.na(x),1, any) else is.na(x)
 if (!any(sample)) warning("data is empty after triming NAs.")
 if (start.) s <-min(time(x)[sample])
 else       s <-start(x)
 if (end.)   e <-max(time(x)[sample])
 else       e <-end(x)
 tfwindow(x,start=s, end=e, warn=F)
}


###############################################

# Non-time dimension methods

###############################################



nseries <- function(x) {UseMethod("nseries")} 
nseries.default <- function(x)  {if (is.matrix(x)) ncol(x) else 1} 

   

 series.names <- function(x)       UseMethod("series.names")
"series.names<-" <- function(x, value)UseMethod("series.names<-")

 series.names.default <- function(x)
   {if (is.null(x)) return(NULL)
    names <- attr(x, "series.names")
    if (is.null(names)) names <- dimnames(x)[[2]]
    if (is.null(names)) names <- paste("Series", seq(ncol(x)))
    names
   }

"series.names<-.default" <- function(x, value){attr(x,"series.names")<-value; x}



select.series <- function(x, series=seqN(nseries(x))) UseMethod("select.series")

select.series.default <- function(x, series=seqN(nseries(x))) {
  names <- series.names(x)
  if (is.character(series)) series <- match(names,series, nomatch=0) > 0
  if(all(0==series) | is.null(series)) r <- NULL
  else {
    r <- classed(tframed(x[, series, drop = F], tframe(x)), class(x))# reconstructor
    series.names(r) <- names[series]
    }
  r
  }

# possibly there should be more attempt to preserve attributes in 
#  select.series but there are problems?:
#     at <- attributes(x)
#     atn <- names(at)
#     atl <- (atn != "dim") & (atn != "dimnames")
#     atn <- atn[atl]
#     at[[!atl]] <- NULL
#     r <- x[,series,drop=F] 
#     list.add(attributes(r), atn) <- at
     

 

tbind <- function(x, ..., pad.start=T, pad.end=T, warn=T)  {UseMethod("tbind")}

tbind.default <- function(x, ..., pad.start=T, pad.end=T, warn=T)
 {# this should work for old tsp vectors and matrices
  if (is.null(x)) stop("first argument cannot be NULL.")
  fr <- frequency(x)
  for (i in list(...)) {if (!is.null(i) && (fr != frequency(i)))
     stop("frequencies must be the same.")}
  fr <- c(fr,1)
  st <- fr %*% start(x) 
  for (i in list(...)) if (!is.null(i)) st <- min(st, fr %*% start(i) )
  en <- fr %*% end(x)
  for (i in list(...)) if (!is.null(i)) en <- max(en, fr %*% end(i) )
  r <- NULL
  # series names (sn) and names/dimnames (nm) do the same thing and sometimes
  # conflict. It is tempting to eliminate nm here, but ...
  sn <- NULL
  nm <- attr(x, "names")
  attr(x, "names") <- NULL
  for (z in append(list(x),list(...)))
   {if (!is.null(z))
    {if (is.matrix(z))
       {if (st == (fr %*% start(z))) before <- NULL
        else  before <-matrix(NA, (fr %*% start(z))-st, dim(z)[2])     
        if (en == (fr %*% end(z))) aft <- NULL
        else  aft    <-matrix(NA, en - (fr %*% end(z)), dim(z)[2])
        r <- cbind(r, rbind( before, z, aft) )
       }
     else 
       {if (st == (fr %*% start(z))) before <- NULL
        else  before <-rep(NA, (fr %*% start(z))-st)     
        if (en == (fr %*% end(z))) aft <- NULL
        else  aft <- rep(NA, en - (fr %*% end(z)))
        r <- cbind(r, c( before, z, aft) )
       }
     sn <- c(sn,series.names(z))
   }}
  if (!is.null(nm)) dimnames(r) <- list(nm,NULL)
  if (length(sn) == ncol(r)) series.names(r) <- sn
  r <- tframed(r, list(start=c((st-1)%/%fr[1], 1+(st-1)%%fr[1]), 
                       frequency=fr[1]))
  if (!(pad.start & pad.end)) r <- trim.na(r, start.=!pad.start, end.=!pad.end)
  if (is.null(r)) warning("intersection is NULL")
  r
 }

############################################################################

#   miscellaneous time calculations  <<<<<<<<<<
#   (Useful utilities not strictly part of tframe)

############################################################################

add.date <- function(date, periods, freq)
  {if (is.null(periods)) periods <- 0
   c(date[1]+(date[2]+periods-1)%/%freq, 1+(date[2]+periods-1)%%freq)
  }


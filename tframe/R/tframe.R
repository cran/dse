
classed <- function(x, cls) {class(x) <- cls; x}

###########################################################################

# freeze and its default method really belong with tfpadi (dsepadi) but the
#  generic and default are needed more generally than the database interface,
#  so they are included here.

freeze <- function(data, ...) UseMethod("freeze") 
 
freeze.default <- function(data, ...){
#  (... further arguments, currently disregarded)
 if ("character"==mode(data)) freeze(tfPADIdata(data, server="ets")) else data} 

#internal utility
# Use this with "for (i in seq(length=m) )" as m==0 returns NULL and for does no loops
seqN <- function(N) {if (0==length(N)) NULL else if (N<=0) NULL else seq(N)}



# start, end, frequency, time need to be masked from R base so that
# tframe methods can work on the tframe attribute rather than class(x)

# The .tframe methods are "default" methods for tframes. Other more specific
#  methods can be defined (see eg start.tstframe for tframes from ts objects). 

# the functions start, end, and frequency in tframe and dse do not 
#  need "...", but the generic in R has it, so it is added here.

start <- function(x, ...) 
 if(is.Ttframed(x)) start(tframe(x), ...) else UseMethod("start") 

start.tframe <- function(x, ...) c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))
#  (... further arguments, currently disregarded)



end <- function(x, ...) 
 if(is.Ttframed(x)) end(tframe(x), ...) else UseMethod("end") 

end.tframe <- function(x, ...)
 c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))
#  (... further arguments, currently disregarded)


frequency <- function(x, ...) 
 if(is.Ttframed(x)) frequency(tframe(x), ...) else UseMethod("frequency") 

frequency.tframe <- function(x, ...) x[3]
#  (... further arguments, currently disregarded)


time <- function(x, ...) 
 if(is.Ttframed(x)) time(tframe(x,...)) else UseMethod("time") 

time.tframe <- function(x, ...) tframed(x[1]+(seq(periods(x))-1)/x[3], tf=x)
#  (... further arguments, currently disregarded)

# periods should give the number of data points in the time direction.
# for consistency check periods needs to look at the data not the tframe,
# i.e. the number of (vector) observations.

periods <- function(x) UseMethod("periods")

periods.default <- function(x) {if (is.array(x)) dim(x)[1] else length(x) }

periods.tframe <- function(x) {1+round((x[2]-x[1])*x[3])}



diff <- function(x, ...) 
 {if (is.Ttframed(x)) {
    tf <- diff(tframe(x),...) 
    tframe(x) <- NULL
    x <- diff(x, ...) 
    return(tframed(x, tf))
    }
  else UseMethod("diff") }

diff.tframe <- function (x, lag = 1, differences = 1, ...) 
 {d <- lag * differences
  tfTruncate(x, start=if(d >= 0) 1+d else NULL, 
                  end=if(d <  0) periods(x)-d else NULL)
 }

#  tfplot and tfprint below provide generic methods for ploting and printing
#  tf time series objects. Plot methods will probably to some processing
#  and eventually call tfplot.default.

tfplot <- function(x, ...)  UseMethod("tfplot")

tfplot.default <- function(x, ..., start.=NULL, end.=NULL, tf=NULL,
       series=seq(nseries(x)), Title=NULL, xlab=NULL, ylab=seriesNames(x), 
       graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
 {#  ... before other args means abbreviations do not work, but otherwise
  # positional matching seems to kick in and an object to be plotted gets used
  #  for start..
  if (!is.tframed(x)) UseMethod("plot")
  else
    {if (is.null(start.) && !is.null(tf)) start. <- start(tf)
     if (is.null(end.)   && !is.null(tf)) end.   <- end(tf)
     others <- list(...)
     tfspan <- x
     # this is a kludge to pass the expansion of shorter series to tbind
     #  and then get the overall time span from the resulting tframe.
     if (0 != length(others)) for (d in others) tfspan <- tbind(tfspan , d)
     tfspan <- tframe(tfspan)
     names <- seriesNames(x)
     Ngraphs <- min(length(series), graphs.per.page)
     if(reset.screen) {
        old.par <- par(mfcol = c(Ngraphs, 1), mar=mar)  
        on.exit(par(old.par)) }
#     tf <- tframe(tfwindow(x, start.=start., end.=end.))
# would be nice if this could expand tf (tfwindow only truncates - need a
# replacement that expands too.)
     tf <- tfwindow(tfspan, start. = start., end. = end., warn = FALSE)
     if (! is.tframe(tf)) browser()
     for (i in series)
       {z <-  matrix(NA, periods(tf), 1+length(others))
        j <- 0
        for (d in append(list(x),others))
    	  {d <- tfwindow(d, tf=tf, warn=FALSE)
	   if (!is.matrix(d)) d <- tframed(as.matrix(d), tf)
	   if(mode(i)=="character") i <- match(i, names)
	   j <- j + 1
	   #z[,j] <- select.series(d, series=i)
	   # tbind in the next will pad d if necessary  to fit tf
	   z[,j] <- select.series(tbind(d,time(tf)), series=i)
	  }
	 tfOnePlot(tframed(z, tf), xlab=xlab, ylab=ylab[i])
	}
    }
  invisible()
 }
 
tfOnePlot <- function(x, xlab=NULL, ylab=NULL, ...)
 {if (!is.tframed(x)) UseMethod("plot")
  else
    {tline <- time(x)
     if(is.null(xlab)) xlab <- ""
     if(is.null(ylab)) ylab <- paste(seriesNames(x), collapse="  ")
     matplot(tline, x, type="l", xlab=xlab, ylab=ylab, ...)
    }
  invisible()
 }

# Note tfprint prints the data. tframePrint  prints the tframe info. 

tfprint <- function(x, ...)  UseMethod("tfprint")

tfprint.default <- function(x,...)
 {xx <- x
  if(1 == nseries(xx)) xx <- matrix(xx, length(xx), 1)
  dimnames(xx) <- list(format(time(tframe(x))), seriesNames(x))
  tframe(xx) <- NULL
  seriesNames(xx) <- NULL
  print(xx, ...)
  invisible(x)
 }



tfwindow <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
  UseMethod("tfwindow")

tfwindow.default <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
  {# With the default warn=TRUE warnings will be issued if no truncation takes
   #  place because start or end is outside the range of data.
   if (is.null(start.) && !is.null(tf)) start. <- start(tf)
   if (is.null(end.)   && !is.null(tf)) end.   <- end(tf)
   # kludge
   x <- ts(x, start=start(x), end=end(x), frequency=frequency(x))
   if (!warn) 
     {opts <- options(warn = -1)
      on.exit(options(opts))
     }
   y <- window(x, start=start., end=end.)
   if (is.matrix(x) && !is.matrix(y) ) y <- matrix(y, length(y), ncol(x))
   y <- tframed(unclass(y), tframe(y))
   seriesNames(y) <- seriesNames(x)
   y
  }


# window a tframe
tfwindow.tframe <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
      tframe(tfwindow(time(x), start.=start., end.=end., tf=tf, warn=warn))

###############################################

#  tframe  methods   <<<<<<<<<<<<

################################################
is.tframe <- function(tf) inherits(tf, "tframe")
is.tframed <- function(x) inherits(tframe(x), "tframe")
# above does not distinguish "true" tframed objects since tframe(x) needs
# to try very hard to return tframes from ts and old tsp objects.
is.Ttframed <- function(x) {!is.null(attr(x, "tframe"))}


tframe <- function(x) UseMethod("tframe")

tframe.default <- function(x){ #extract the tframe
	if(is.null(x)) NULL 
	else if (!is.null(attr(x, "tframe"))) attr(x, "tframe") # extractor
	else if (!is.null(tsp(x)))    classed(tsp(x), "tframe") # extractor
	else if(is.vector(x)) classed(c(1,length(x),1), "tframe") # extractor
	else if(is.matrix(x)) classed(c(1,  nrow(x),1), "tframe") # extractor
	else if(is.array(x) ) classed(c(1,dim(x)[1],1), "tframe") # extractor
	else NULL
}

# using classed(tsp(as.ts(x)), "tframe") in the last line above 
# makes too many things into tframes (eg lists)

"tframe<-" <- function(x, value) UseMethod("tframe<-")

"tframe<-.default" <- function(x, value) {
  # It is tempting in the next to try and make a ts if value is from a ts, 
  #  but that will not work for cases were x does not fit the ts model, so
  #  that would break  tframe(x) <- tframe(y) 
  if(is.null(value)) tf <- value 
  else if(is.tframe(value)) tf <- value 
  else if (is.numeric(value) && length(value) == 3) tf <- value
  else if (is.list(value))
 	{if( is.null(value$start) & is.null(value$end) )
 	    stop("value must be a tframe or a list of arguments for ts().")
 	 attr(x, "tframe") <- NULL #Splus 3.3 bug work around
	 tf <- tframe(do.call("ts", append(list(x), value)))
 	}
  # next is checking after the fact, but value may just be start and freq
  #  which is not enough to know periods
  attr(x, "tframe") <- tf
  if((!is.null(value)) && !check.tframeConsistent(tframe(x), x))
     stop("time frame value in tframe assignment is not consistent with data.")
  x
}

# making tframed generic allows tframed.TSdata to specify input and output names
# but there may be a better way and leave "tframe<-" as the real generic.

tframed <- function(x, tf=NULL, names = NULL) UseMethod("tframed")

tframed.default <- function(x, tf=NULL, names = NULL){
  # return x as a tframed object with tframe tf
  if (!is.null(names))  seriesNames(x) <-  names
  if (is.null(tf)) tf <- tframe(x) # this generates a default
  tframe(x) <- tf
  x
 }


###############################################

#  Generic .tframe methods (these act on the tframe not on the data)

###############################################


tfprint.tframe <- function(x, ...) UseMethod("tframePrint")
tframePrint <- function(x, ...) UseMethod("tframePrint")

tframePrint.default <- function(x, digits=NULL, quote=TRUE, prefix="", ...) 
  {if (! is.tframe(x)) x <- tframe(x)
   invisible(print(unclass(x), quote=quote, prefix=prefix, ...)) }



tfTruncate.tframe <- function(x, start=NULL, end=NULL)
    {# like window but uses indexes rather than dates 
     if (!is.null(end))   x[2] <- x[1] + (end-1)/x[3]
     if (!is.null(start)) x[1] <- x[1] + (start-1)/x[3]
     x
    }



tfExpand.tframe <- function(x, add.start=0, add.end=0)
    {x[2] <- x[2] + add.end/x[3]
     x[1] <- x[1] - add.start/x[3]
     x
    }


check.tframeConsistent <- function(tf, x) UseMethod("check.tframeConsistent")

check.tframeConsistent.default <- function(tf, x) periods(tf) == periods(x)

test.Equaltframes <- function(tf1, tf2) UseMethod("test.Equaltframes")

test.Equaltframes.default <- function(tf1, tf2) { all(tf1==tf2)}



# Following could be used to do date comparisons like start() < end()


earliestStart <- function(x, ...)
    start(append(list(x),list(...))[[earliestStartIndex(x, ...)]])

earliestStartIndex <- function(x, ...) UseMethod("earliestStartIndex")

earliestStartIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("earliestStartIndex", tf) #dispatch on 1st element of tf
  }

earliestStartIndex.tframe <- function(x, ...) 
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




earliestEnd <- function(x, ...)
    end(append(list(x),list(...))[[earliestEndIndex(x, ...)]])

earliestEndIndex <- function(x, ...) UseMethod("earliestEndIndex")

earliestEndIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("earliestEndIndex", tf) #dispatch on 1st element of tf
  }

earliestEndIndex.tframe <- function(x, ...) 
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



latestStart <- function(x, ...)
    start(append(list(x),list(...))[[latestStartIndex(x, ...)]])

latestStartIndex <- function(x, ...) UseMethod("latestStartIndex")

latestStartIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("latestStartIndex", tf)
  }


latestStartIndex.tframe <- function(x, ...) 
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



latestEnd <- function(x, ...)
    end(append(list(x),list(...))[[latestEndIndex(x, ...)]])

latestEndIndex <- function(x, ...) UseMethod("latestEndIndex")

latestEndIndex.default <- function(x, ...)
  {tf <- list(tframe(x))
   for (i in list(...)) tf <- append(tf, list(tframe(i)))
   do.call("latestEndIndex", tf)
  }

latestEndIndex.tframe <- function(x, ...) 
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

#   .rts cts and its methods

###############################################

"tframe<-.rts" <- function(x, value) {rts(x) <- value; x}
"tframe<-.cts" <- function(x, value) {cts(x) <- value; x}
"tframe<-.its" <- function(x, value) {its(x) <- value; x}


###############################################

#  stamped specific methods   <<<<<<<<<<<<
#  stamped class TS have a date/time stamp associated with each time point
################################################

#check.tframeConsistent.stamped <- function(tf, x)
#  {periods(x) == periods(tf)}

test.Equaltframes.stamped <- function(tf1, tf2)
  {all(tf1$stamp == tf2$stamp)}

periods.stamped <- function(x) length(tframe(x))

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

if (!exists("lag")) lag <- function(x, ...) UseMethod("lag")

if (!exists("lag.default"))  lag.default <- function(x, ...) {stop("no lag function") }




splice <- function(mat1, mat2, ...) UseMethod("splice")

splice.default <- function(mat1, mat2, ...)
{#  (... further arguments, currently disregarded)
 # splice together 2 time series matrices. If data  is provided in both for
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
 r1 <- tframed(r1, list(start=earliestStart(mat1,mat2), 
                        end =latestEnd(mat1,mat2), frequency=freq))
 r1
}


#if( !exists("tsmatrix.default"))  
#  {if(exists("tsmatrix")) tsmatrix.default <- tsmatrix 
#   else tsmatrix.default <- function(x, ...) 
#	     {tbind(x, ..., pad.start=FALSE, pad.end=FALSE) }
#  }
#
#tsmatrix <- function(x, ...)
# {# the default tsmatrix messes up because it gets some time info. (from
#  #  start or end) but not tsp info.
#  if (is.Ttframed(x)) tbind(x, ..., pad.start=FALSE, pad.end=FALSE)
#  else 
#    {#warning("Using tsmatrix which should be defunct. Consider using tbind and tframe methods.")       
#     tsmatrix.default(x,  ...)
#    }
# }


tfTruncate <- function(x, start=NULL, end=NULL) UseMethod("tfTruncate")
  # similar to window but start and end specify periods relative to the 
  #   beginning (eg x[start:end] for a vector).
  #   NULL means no truncation.


tfTruncate.default <- function(x, start=NULL, end=NULL)
    {tf <- tfTruncate(tframe(x), start, end)
     if (is.null(start)) start <- 1
     if (is.matrix(x)) 
        {if (is.null(end)) end <- dim(x)[1]
         z <- x[start:end,,drop=FALSE]
        }
     else 
        {if (is.null(end)) end <- length(x)
         z <- x[start:end]
        }
     tframe(z) <- tf
     z
    }

tfExpand <- function(x, add.start=0, add.end=0) UseMethod("tfExpand")
  # expand (a tframe) by add.start periods on the beginning
  # and add.end periods on the end

tfExpand.default <- function(x, add.start=0, add.end=0)
    {tf <- tfExpand(tframe(x), add.start=add.start, add.end=add.end)
     select.series(tbind(x,time(tf)), series=1)
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
 s <- if (start.) min(time(x)[sample]) else start(x)
 e <- if (end.)   max(time(x)[sample]) else end(x)
 tfwindow(x, start.=s, end.=e, warn=FALSE)
}


###############################################

# Non-time dimension methods

###############################################



nseries <- function(x) UseMethod("nseries") 
nseries.default <- function(x)  {if (is.matrix(x)) ncol(x) else 1} 

   

 seriesNames <- function(x)       UseMethod("seriesNames")
"seriesNames<-" <- function(x, value)UseMethod("seriesNames<-")

 seriesNames.default <- function(x)
   {if (is.null(x)) return(NULL)
    names <- attr(x, "seriesNames")
    if (is.null(names)) names <- dimnames(x)[[2]]
    if (is.null(names)) names <- paste("Series", seq(ncol(x)))
    names
   }

"seriesNames<-.default" <- function(x, value){attr(x,"seriesNames")<-value; x}



select.series <- function(x, series=seqN(nseries(x))) UseMethod("select.series")

select.series.default <- function(x, series=seqN(nseries(x))) {
  names <- seriesNames(x)
  if (is.character(series)) series <- match(names,series, nomatch=0) > 0
  if(all(0==series) | is.null(series)) r <- NULL
  else {
#    r <- classed(tframed(x[, series, drop = FALSE], tframe(x)), class(x))# reconstructor
#   tframe assignment cannot guarantee that the object has the right structure
#   for a class, so above can give a deformed object in the class.
    r <- tframed(x[, series, drop = FALSE], tframe(x))
    seriesNames(r) <- names[series]
    }
  r
  }


tbind <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
   UseMethod("tbind")

tbind.default <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
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
     sn <- c(sn,seriesNames(z))
   }}
  if (!is.null(nm)) dimnames(r) <- list(nm,NULL)
  if (length(sn) == ncol(r)) seriesNames(r) <- sn
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


############################################################

#  Utility function for time series noise   <<<<<<<<<<

############################################################

makeTSnoise <- function(sampleT,p,lags,noise=NULL, rng=NULL,
                        SIGMA=NULL, sd=1, noise.model=NULL,
                        noise.baseline=0,
                        tf=NULL, start=NULL,frequency=NULL)
 {# CAUTION: changes here can affect historical comparisons.
  # noise.baseline is added to noise. It should be either a scalar, a matrix of
  #   the same dimension as noise (or noise generated by noise.model), or a
  #   vector of length equal to the dimension of the noise process (which will
  #   be replicated for all periods.)
 require("setRNG")
 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }

  if ( (!is.null(noise)) & (!is.null(noise.model)) )
    stop("noise and noise.model cannot both be specified.")

  if(!is.null(noise.baseline) && is.matrix(noise.baseline) &&
    (sampleT < dim(noise.baseline)[1]))
      {warning("sampleT (and start date) for noise adjusted to match noise.baseline")
       sampleT <- dim(noise.baseline)[1]
      }

 # Note: noise is added to initial conditions.
 if (!is.null(noise.model))
   {require("dse1")
    noise <- output.data(simulate(noise.model, sampleT=sampleT+lags))
    noise <- list(w0=noise[1:lags,,drop=FALSE], w=noise[lags+seq(sampleT),,drop=FALSE])
   }
  if (is.null(noise))
   {w0 <-matrix(NA,lags,p)
    w <- matrix(NA,sampleT,p)
    if (!is.null(SIGMA))
       {if(length(SIGMA) == 1) SIGMA <- diag(SIGMA, p)
	W <- t(chol(SIGMA))
        w <- t(W %*% t(matrix(rnorm((lags+sampleT)*p),(lags+sampleT),p)))
	w0 <- w[1:lags,]
	w  <- w[-c(1:lags),]
       }
    else
      {if (length(sd)==1) sd <-rep(sd,p)
       for (i in 1:p)
         {w0[,i] <- rnorm(lags,sd=sd[i])
          w[,i]  <- rnorm(sampleT,sd=sd[i])
         }
      }
    noise <- list(w=w, w0=w0)
   }

  if(!is.null(noise.baseline))
     {if (is.vector(noise.baseline))
        {if(length(noise.baseline)==1) noise$w <- noise$w + noise.baseline
         else if(length(noise.baseline)==1)
           noise$w <-noise$w + t(array(noise.baseline, rev(dim(noise$w))))
         else stop("noise.baseline vector is not correct length.")
        }
      else  noise$w <- noise$w + noise.baseline
     }

  if(!is.null(tf)) tframe(noise$w) <- tf
  else if(!is.null(start))
      {if (is.null(frequency))
         {frequency <- 1
          warning("start set but frequency not specified. Using frequency=1.")
         }
       else noise$w <-tframed(noise$w, list(start=start, frequency=frequency))
       if (is.tframed(noise.baseline) &&
           test.equal(tframe(noise.baseline),tframe(noise$w)))
           {warning("tframe of noise set to tframe of noise.baseline.")
            tframe(noise$w )<-tframe(noise.baseline)
            if(!all(dimnames(noise$w)[[2]] == dimnames(noise.baseline)[[2]]))
              warning("noise names and noise.baseline names do not correspond.")
           }
      }
  append(noise, list(sampleT=sampleT, rng=rng,
     SIGMA=SIGMA, sd=sd, noise.model=noise.model,version=as.vector(version)))
 }


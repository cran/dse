###############################################

#   .ts methods (for the object) and .tstframe methods (for the tframe)

###############################################

"seriesNames<-.ts" <- function (x, value) 
  {attr(x, "seriesNames") <- value
   if (is.matrix(x)) dimnames(x) <- list(NULL, value)
   x
  }

tframe.ts <- function(x){classed(tsp(x), c("tstframe", "tframe"))} # extractor

"tframe<-.ts" <- function(x, value) {do.call("ts", append(list(x), value))}

select.series.ts <- function(x, series=seqN(nseries(x))) {
  names <- seriesNames(x)
  if (is.character(series)) series <- match(names,series, nomatch=0) > 0
  if(all(0==series) | is.null(series)) r <- NULL
  else {
    r <- x[, series, drop = FALSE]
    seriesNames(r) <- names[series]
    }
  r
  }

tbind.ts <- function(x, ..., pad.start=TRUE, pad.end=TRUE, warn=TRUE)
 {# this is used like old tsmatrix should produce a column matrix from a
  #  single vector
  nm <- seriesNames(x)
  for (z in list(...)) {
    if (!is.null(z)) {
      nm <- c(nm, seriesNames(z))
      if (!is.ts(z)) z <- ts(z,start=start(z),end=end(z),frequency=frequency(z))
      x <- cbind(x, z)
      }
    }
  if (!is.matrix(x)) x <- ts(matrix(x, length(x),1),
                      start=start(x), end=end(x), frequency=frequency(x))
  if (!pad.start | !pad.end)
     x <- trim.na(x, start. = !pad.start, end. = !pad.end)
  seriesNames(x) <- nm
  x
 }

tfwindow.ts <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
  {# With the default warn=T warnings will be issued if no truncation takes
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
   seriesNames(y) <- seriesNames(x)
   y
  }


# The .tstframe methods should work on tframe(ts(whatever))  from ts objects. 

start.tstframe <- function(x) c(floor(x[1]+getOption("ts.eps")),
    round(1 + ((x[1]+getOption("ts.eps"))%%1)*x[3]))
end.tstframe   <- function(x) c(floor(x[2]+getOption("ts.eps")),
    round(1 + ((x[2]+getOption("ts.eps"))%%1)*x[3]))
periods.tstframe <- function(x)  {1+round((x[2]-x[1])*x[3])}
frequency.tstframe <- function(x) x[3]
time.tstframe <- function(x) {x[1] + (seq(periods(x))-1)/x[3]}

tfTruncate.tstframe <- function(x, start=NULL, end=NULL) 
    {if (!is.null(end))   x[2] <- x[1] + (end-1)/x[3]
     if (!is.null(start)) x[1] <- x[1] + (start-1)/x[3]
     x
    }

tfExpand.tstframe <- function(x, add.start=0, add.end=0) 
    {x[2] <- x[2] + add.end/x[3]
     x[1] <- x[1] - add.start/x[3]
     x
    }


earliestStartIndex.tstframe <- function(x, ...) 
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

earliestEndIndex.tstframe <- function(x, ...) 
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

latestStartIndex.tstframe <- function(x, ...) 
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

latestEndIndex.tstframe <- function(x, ...) 
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

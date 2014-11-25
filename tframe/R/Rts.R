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

selectSeries.ts <- function(x, series=seqN(nseries(x))) {
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
     x <- trimNA(x, start. = !pad.start, end. = !pad.end)
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


# The .tstframe methods work on tframe the tframe of a ts object. 
# The .tframe (next) methods should work for most <tstframe tframe> tframes.
# Following are a couple that are slightly different.

start.tstframe <- function(x, ...) c(floor(x[1]+getOption("ts.eps")),
    round(1 + ((x[1]+getOption("ts.eps"))%%1)*x[3]))
#  (... further arguments, currently disregarded)

end.tstframe   <- function(x, ...) c(floor(x[2]+getOption("ts.eps")),
    round(1 + ((x[2]+getOption("ts.eps"))%%1)*x[3]))
#  (... further arguments, currently disregarded)


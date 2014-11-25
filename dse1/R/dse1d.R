#   2000/04/20 14:50:54 
# For installation instructions see the file read.me or the brief user's
#    guide (postscipt file guide.ps).


############################################################

#     Methods which are generic for models and TSdata    <<<<<<<<<<

# generic methods for TS structures (ie with input and output) <<<<<<<<<<
#                   see also tframe.s                          <<<<<<<<<<

############################################################################
# periods and periods.default are defined in tfame.s

input.periods <- function(x, ...)UseMethod( "input.periods")
output.periods <- function(x, ...)UseMethod("output.periods")

input.start <- function(x)UseMethod("input.start")
output.start <- function(x)UseMethod("output.start")

input.end <- function(x)UseMethod("input.end")
output.end <- function(x)UseMethod("output.end")

input.frequency <- function(x)UseMethod("input.frequency")
output.frequency <- function(x)UseMethod("output.frequency")




input.data <- function(x, ...)UseMethod("input.data")
output.data <- function(x, ...)UseMethod("output.data")

input.data.default <- function(x) x$input
output.data.default <- function(x) x$output

"input.data<-" <- function(x,  value)  UseMethod("input.data<-")
"output.data<-" <- function(x,  value) UseMethod("output.data<-")

"input.data<-.default" <- function(x,  value) {x$input <- value; x}
"output.data<-.default" <- function(x,  value){x$output <- value; x}

# The logic (revised as of May26, 1998) is that series.names can be an attribute
# of any object (like a matrix) and has been moved to tframe.



 input.series.names <- function(data)UseMethod( "input.series.names")
output.series.names <- function(data)UseMethod("output.series.names")

 "input.series.names<-" <- function(x, value)UseMethod( "input.series.names<-")
"output.series.names<-" <- function(x, value)UseMethod("output.series.names<-")


############################################################################

#    methods for TSmodels  <<<<<<<<<<
# and check elsewhere too

############################################################################

output.series.names.TSmodel <- function(model)
 {# return output names if available in the object,
  # otherwise return "out" pasted with integers.
  if (!is.null(attr(model, "output.series.names")))
                          return(attr(model, "output.series.names"))
  else if (0 != output.dimension(model)) 
                   return(paste("out", seq(output.dimension(model)), sep=""))  
  else return(NULL)
 }

input.series.names.TSmodel <- function(model)
 {# return input names if available in the object,
  # otherwise return "in" pasted with integers.
  if (!is.null(attr(model, "input.series.names")))
                          return(attr(model, "input.series.names"))
  else if (0 != input.dimension(model)) 
                   return(paste("in", seq(input.dimension(model)), sep=""))  
  else return(NULL)
 }

"series.names<-.TSmodel" <- function(x, value)
   { input.series.names(x) <-  value$input
    output.series.names(x) <-  value$output
    x
   }

"input.series.names<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != input.dimension(x))
       stop("model input dimension and number of names do not match.")
    attr(x, "input.series.names")  <- value
    x
   }

"output.series.names<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != output.dimension(x))
       stop("model output dimension and number of names do not match.")
    attr(x, "output.series.names")  <- value
    x
   }


############################################################################

#    methods for TSestModels  <<<<<<<<<<
# and check elsewhere too   (esp. for start, end and frequency)

############################################################################
periods.TSestModel <- function(x)periods(x$data)
input.periods.TSestModel <- function(x)input.periods(x$data)

input.data.TSestModel <- function(x, ...) input.data(x$data, ...)
output.data.TSestModel <- function(x, ...)output.data(x$data, ...)

input.dimension.TSestModel <- function(x)  input.dimension(x$data)
output.dimension.TSestModel <- function(x) output.dimension(x$data)

input.series.names.TSestModel <- function(obj)
 {m <- input.series.names(obj$model)
  d <- input.series.names(obj$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

output.series.names.TSestModel <- function(obj)
 {m <- output.series.names(obj$model)
  d <- output.series.names(obj$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

series.names.TSestModel <- function(obj)
  {list(input=input.series.names(obj), output=output.series.names(obj))}

"series.names<-.TSestModel" <- function(x, value)
   { input.series.names(x) <- value$input
    output.series.names(x) <- value$output
    x
   }


"input.series.names<-.TSestModel" <- function(x, value)
   {input.series.names(x$model) <- value;
    input.series.names(x$data ) <- value;
    x
   }
"output.series.names<-.TSestModel" <- function(x, value) 
   {output.series.names(x$model) <- value;
    output.series.names(x$data ) <- value;
    x
   }

identifiers.TSestModel <- function(obj){identifiers(TSdata(obj))}
sourcedb.TSestModel <- function(obj){sourcedb(TSdata(obj))}
sourceserver.TSestModel <- function(obj){sourceserver(TSdata(obj))}
source.info.TSestModel <- function(obj){source.info(TSdata(obj))}


############################################################################

#    methods for TSdata  <<<<<<<<<<
#      note there are some methods in dse1c too

############################################################################

is.TSdata <- function(obj) { inherits(obj, "TSdata")}

print.TSdata <- function(x, ...)
{  if(0 != (input.dimension(x)))
     {cat("input data:\n")
      print(input.data(x),...)
      if(!is.null(x$input.transformations))
          {cat("input.transformations:\n")
           print(x$input.transformations, ...)
          }
      if(!is.null(input.series.names(x)))
          {cat("input.names:\n")
           print(input.series.names(x), ...)
          }
      cat("\n")
     }
  if(0 != (output.dimension(x)))
     {cat("output data:\n")
      print(output.data(x),...)
      if(!is.null(x$output.transformations))
         {cat("output.transformations:\n")
          print(x$output.transformations, ...)
         }
      if(!is.null(output.series.names(x)))
         {cat("output.names:\n")
          print(output.series.names(x), ...)
         }
     }
   cat("\n")
   if(!is.null(x$retrieval.date))
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) 
     {cat("source:\n")
      print(x$source)
     }
  invisible(x)
}

summary.TSdata <- function(object, ...)
  {d <- output.data(object)
   if (!is.tframed(d)) d <- as.ts(d)
   st <- start(d)
   en <- end(d)
   fr <- frequency(d)
   d <- cbind(d,input.data(object))

   classed(list(  # summary.TSdata constructor
      description=object$description,
      start=st,
      end=en,
      freq=fr,
      sampleT= nrow(output.data(object)),
      p=output.dimension(object),
      m=input.dimension(object),
      ave=apply(d,2,mean),
      max=apply(d,2,max),
      min=apply(d,2,min),
      retrieval.date=object$retrieval.date, 
      source=object$source), "summary.TSdata")
  }

print.summary.TSdata <- function(x, digits=options()$digits)
{  if (!is.null(x$description)) cat(x$description,"\n")
   cat("start.: ", x$start, " end.: ",  x$end," Frequency: ", x$freq,"\n")
   cat("sample length = ",x$sampleT,"\n")
   cat("number of outputs=",x$p, "   number of inputs=",x$m, "\n")
   cat("average :\n")
   print(x$ave)
   cat("max    :\n")
   print(x$max)
   cat("min    :\n")
   print(x$min)
   cat("\n")
   if(!is.null(x$retrieval.date)) 
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) {cat("source:\n"); print(x$source)  }
   invisible(x)
}


tfplot.TSdata <- function(..., start.=NULL,end.=NULL, Title="", reset.screen=T,
        select.inputs =seq(length= input.dimension(data)),
        select.outputs=seq(length=output.dimension(data)),
	mar=if(is.R()) c(3.1,6.1,1.1,2.1) else c(5.1,6.1,4.1,2.1) ,
        graphs.per.page=5, ylab=NULL)
{# plot input and output data.
 # ... is a list of objects of class TSdata (with similar input
 #       and output dimensions.
 # Note that using ... like this means it cannot be used to pass additional
 #   arguments to plot, so unfortunately all necessary plot arguments must be 
 #   explicit in the arguments to tfplot.TSdata
 # start. is the starting point (date)  and end. the ending point for
 # plotting. If not specified the whole sample is plotted.
 # use dev.ask(T) to pause before a new page is printed
  data <-list(...)[[1]]
  if (!is.TSdata(data)) TS.error.exit(clss="TSdata") 
  data <- freeze(data)
  Ngraphs <- length(select.outputs) + length(select.inputs)
  if(reset.screen)
    {Ngraphs <- min(Ngraphs, graphs.per.page)
     old.par <- par(mfcol = c(Ngraphs, 1), mar=mar) # previously c(5.1,6.1,4.1,2.1)) 
     on.exit(par(old.par))
    }
  if (!is.null(ylab))
       names <- rep(ylab,input.dimension(data) + output.dimension(data))
  else names <-  c(   input.series.names(data),  output.series.names(data))
#  if (0 != length(select.inputs)) 
    {for (i in select.inputs) 
      {j <- 0
       z <- matrix(NA, input.periods(data), length(list(...)))
       if(mode(i)=="character") i <- match(i, input.series.names(data))
       for (data in list(...) ) 
         {if (!is.TSdata(data))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
          data <- freeze(data)
#         Rbug 0.65 in ts   z <- cbind(z, input.data(data, series = i))
          j <- j + 1
          z[, j] <- input.data(data, series = i)
         }
       tframe(z) <-tframe(input.data(data))
       if (!is.null(start.)) z <- tfwindow(z,start=start.)
       if (!is.null(end.))   z <- tfwindow(z,end=end.)
       tfplot(z,ylab=names[i])  # tsplot
       if(i==1) title(main = Title)
    } }
  for (i in select.outputs) 
    {j <-0
     if(mode(i)=="character") i <- match(i, output.series.names(data))
     z <- matrix(NA,output.periods(data),length(list(...)))
     #z <-NULL 
     for (data in list(...) ) 
       {if (!is.TSdata(data))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
        data <- freeze(data)
        j <- j+1
        z[,j]<-output.data(data, series=i) 
       }
     tframe(z) <-tframe(output.data(data))
     if (!is.null(start.)) z <- tfwindow(z,start=start.)
     if (!is.null(end.))   z <- tfwindow(z,end=end.)
     tfplot(z,ylab=names[input.dimension(data) + i]) # tsplot
     if((0 == input.dimension(data)) & (i==1)) title(main = Title)
    }
  invisible()
}
 

#was at one point to avoid Rbug: tframed(unclass(x$output)[,series,drop=F], tframe(x$output))
#  previously else tframed(x$input[,series,drop=F], tframe(x$input))


input.data.TSdata <- function(x, series=seq(input.dimension(x)))
  {if (is.null(x$input)) NULL  else select.series(x$input, series=series) }

output.data.TSdata <- function(x, series=seq(output.dimension(x)))
  {if (is.null(x$output)) NULL  else select.series(x$output, series=series) }



# Rbug setting the class in the following two functions was 
#  necessary for NULL assignments
# cls <- class(x); x$input <-newinput; class(x) <- cls
# instead of simply
# x$input <-newinput

"input.data<-.TSdata" <- function(x, value) 
   {cls <- dseclass(x); x$input <- value; dseclass(x) <- cls
#    input.series.names(x) <- dimnames(x)[[2]]
    x
   }


"output.data<-.TSdata" <- function(x, value)
   {cls <- dseclass(x); x$output <-value; dseclass(x) <- cls
# not nec. with change to attr    output.series.names(x) <- dimnames(x)[[2]]
    x
   }

# Note: series names changed to an attribute of input and output for 
#     data but not for models!!!!

series.names.TSdata <- function(x)
 {list(input=input.series.names(x), output=output.series.names(x))}

 input.series.names.TSdata <- function(x) {series.names( input.data(x))}
output.series.names.TSdata <- function(x) {series.names(output.data(x))}

"series.names<-.TSdata" <- function(x, value)
   { input.series.names(x) <-  value$input
    output.series.names(x) <-  value$output
    x
   }

"input.series.names<-.TSdata" <- function(x,  value)
   {if (length( value) != input.dimension(x))
       stop("number of series and number of names do not match.")
    attr(input.data(x), "series.names")  <- value
    x
   }
"output.series.names<-.TSdata" <- function(x,  value) 
    {if (length( value) != output.dimension(x))
       stop("number of series and number of names do not match.")
    attr(output.data(x), "series.names")  <- value
    x
   }

input.dimension.TSdata <- function(obj)
   {if (is.null(obj$input)) 0 else nseries(obj$input)}

output.dimension.TSdata <- function(obj)
   {if (is.null(obj$output)) 0 else nseries(obj$output)}

start.TSdata <- function(data)
{i  <-  input.start(data)
 o  <- output.start(data)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 else return(c(i,o))
}

input.start.TSdata <- function(data)
 {if (is.null(data$input))  return(NULL)
  else return(start(data$input))
 }

output.start.TSdata <- function(data)
 {if (is.null(data$output))  return(NULL)
  else return(start(data$output))
 }

end.TSdata <- function(data)
{i  <-  input.end(data)
 o  <- output.end(data)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

input.end.TSdata <- function(data)
 {if (is.null(data$input))  return(NULL)
  else return(end(data$input))
 }

output.end.TSdata <- function(data)
 {if (is.null(data$output))  return(NULL)
  else return(end(data$output))
 }

frequency.TSdata <- function(data)
{i  <-  input.frequency(data)
 o  <- output.frequency(data)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

input.frequency.TSdata <- function(data)
 {if (is.null(data$input))  return(NULL)
  else return(frequency(data$input))
 }

output.frequency.TSdata <- function(data)
 {if (is.null(data$output))  return(NULL)
  else return(frequency(data$output))
 }


periods.TSdata <- function(data) UseMethod("output.periods")
output.periods.TSdata <- function(data)  dim(output.data(data))[1]
input.periods.TSdata <- function(data)  dim(input.data(data))[1]

tbind.TSdata <- function(d1, d2)
 {if( ! (is.TSdata(d1) & is.TSdata(d2)))
     stop("tbind requires arguments to be of a similar type (ie. TSdata).")
  if ((0 != input.dimension(d1)) || (0 != input.dimension(d2)) )
    input.data(d1) <- tbind(input.data(d1),input.data(d2))
  if ((0 != output.dimension(d1)) || (0 != output.dimension(d2)) )
    output.data(d1) <- tbind(output.data(d1),output.data(d2))
  d1
 }



test.equal.TSdata <- function(d1,d2, fuzz=1e-16)
  {r <- T
   if (r & (!is.null(d1$input)))
     {if(is.null(d2$input)) r <- F
      else  r <-test.equal.matrix(d1$input, d2$input, fuzz=fuzz)
     }
   if (r & (!is.null(d1$output)))
     {if(is.null(d2$output)) r <- F
      else r <-test.equal.matrix(d1$output, d2$output, fuzz=fuzz)
     }
   r
  }


identifiers.TSdata <- function(obj) {identifiers(obj$source)}
sourcedb.TSdata <- function(obj) {sourcedb(obj$source)}
sourceserver.TSdata <- function(obj) {sourceserver(obj$source)}
source.info.TSdata <- function(obj) {source.info(obj$source)}




TSdata <- function (data=NULL, ...) {UseMethod("TSdata")} 
 

TSdata.default <- function(data=NULL, input=NULL, output=NULL, ...)  
{if (is.null(data) && (!is.null(input) | !is.null(output) ))
    data <- classed(list(input=input, output=output), "TSdata") # constructor
 else 
    data <- classed(data, "TSdata")   # constructor keeps other list elements
  
 if(!is.list(data)) stop("TSdata.default could not construct a TSdata format.")	    
 if   ( 0 == input.dimension(data)   &    0 == output.dimension(data) 
    |  (0 != output.dimension(data) && !is.matrix(output.data(data)))
    |  (0 !=  input.dimension(data) && !is.matrix( input.data(data))) )
      stop("TSdata.default could not construct a TSdata format.")
 data
}

TSdata.TSdata <- function(data, ...)  {data} # already TSdata 

TSdata.TSestModel <- function(data, ...)  {data$data} # extract TSdata 

as.TSdata <- function(d) 
 {# Use whatever is actually $input and $output,
  #strip any other class, other parts of the list, and DO NOT use
  # input.data(d) and output.data(d) which may eg reconstitute something
  TSdata(input=d$input, output=d$output)
 }

tframed.TSdata <- function(x,  tf=NULL, input.names=NULL, output.names=NULL)  
{# switch to tframe representation
 if(0 != (output.dimension(x)))
       output.data(x) <-tframed(output.data(x), tf=tf, names=output.names)
 if (0 != (input.dimension(x)))
        input.data(x) <-tframed(input.data(x),  tf=tf, names= input.names)
 x
}  


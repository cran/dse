# For installation instructions see the file read.me or the brief user's
#    guide (postscipt file guide.ps).

##############################################################################

############################################################################

#    functions for DSE interface to Time Series Protocol for      <<<<<<<<<<
#      Application Database Interface (TSPADI) data interface    <<<<<<<<<<

############################################################################

# Functions in this file now handle only the TSdata aspects. The database
#  interface has been pushed down into a time series matrix class "tfPADIdata"
#  defined in the file tfpadi.

# It is not clear that the c("TSPADIdata", "TSdata") class is necessary,
#  as the tfPADIdata should be just another time series matrix class,
#  but in an attempt to smooth the transition that is not being attempted 
#  in one step. Also, the name tfPADIdata might better be TSPADIdata,
#  but that would be confusing for the time being.

# Also, freeze is used to do both freeze and TSdata since these have been 
#  done together previously. This should be cleaned up.

############################################################

#   Definition of class c("TSPADIdata", "TSdata") <<<<<<<<<<

############################################################




# This TSPADIdata constructor was once called make.TSnames

TSPADIdata <- function( output=NULL,           input=NULL,
                        output.server=server,  input.server=server,
                        output.db=db,          input.db=db,
                        output.transforms="",  input.transforms="", 
                        output.names=NULL,     input.names=NULL,
                         start=NA, end=NA, frequency=NA, 
                         pad=FALSE, pad.start=pad, pad.end=pad,
                         server="", db="", start.server=NULL, 
                         server.process=NULL, cleanup.script=NULL,
                         stop.on.error=T, warn=T)
  {i <- if (is.null(input)) NULL else tfPADIdata(input, 
      transforms=input.transforms, names=input.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

   o <- if (is.null(output)) NULL else tfPADIdata(output, 
      transforms=output.transforms, names=output.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end,
      server=output.server, db=output.db, start.server=start.server, 
      server.process=server.process, cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

    classed(list(input=i, output=o), c("TSPADIdata", "TSdata")) # constructor 
   }





TSPADIdata2 <- function(input=NULL, output=NULL,
    start = NA, end = NA, frequency = NA, pad.start = FALSE, 
    pad.end = FALSE,  start.server = NULL, 
    server.process = NULL, cleanup.script = NULL, stop.on.error = T, 
    warn = T)
  {i <- o <- NULL
   for (j in seq(length=length( input))) i <- cbind(i,  input[[j]])
   for (j in seq(length=length(output))) o <- cbind(o, output[[j]])
   TSPADIdata(input =     i[3,], output=            o[3,],
        input.server=     i[1,], output.server=     o[1,],
        input.db=         i[2,], output.db=         o[2,],
        input.transforms= i[4,], output.transforms= o[4,],
        input.names=      i[5,], output.names=      o[5,],
      start = start, end = end, frequency = frequency, pad.start = pad.start, 
      pad.end = pad.end, start.server = start.server, 
      server.process = server.process, cleanup.script = cleanup.script,
      stop.on.error = stop.on.error, warn = warn)
   }
  



modify.TSPADIdata <- function(obj,
                        output=NA,             input=NA,
                        output.server=NA,      input.server=NA,
                        output.db=NA,          input.db=NA,
                        output.transforms=NA,  input.transforms=NA, 
                        output.names=NA,       input.names=NA,
                         start=NA, end=NA, frequency=NA, 
                         pad.start=NA, pad.end=NA,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA)
  {
   if( (!all(is.na(c(input, input.server, input.db, input.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    input.data(obj) <-  modify.tfPADIdata(input.data(obj),
      series=input, transforms=input.transforms, names=input.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

   if( (!all(is.na(c(output, output.server, output.db, output.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    output.data(obj) <-  modify.tfPADIdata(output.data(obj),
      series=output, transforms=output.transforms, names=output.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)
    obj
   }



############################################################

#     methods for TSPADIdata class objects <<<<<<<<<<

############################################################

print.TSPADIdata <- function(x, ...) {print.default(x) }

is.TSPADIdata <- function(obj) {inherits(obj, "TSPADIdata") }

# TSdata methods should work for start, end, frequency


tsp.TSPADIdata <- function(x)
  {i <- tsp( input.data(x))
   o <- tsp(output.data(x))
   if (is.null(o)) return(i)
   if (is.null(i)) return(o)
   if (!all(i == o)) 
      warning("tsp results differ for input and output data. Using output")
   o
}

input.periods.TSPADIdata <- function(data) periods( input.data(data))  
output.periods.TSPADIdata <- function(data) periods(output.data(data))  
periods.TSPADIdata <- function(data) periods(output.data(data))

 
input.data.TSPADIdata <- function(x, series=seq(length=input.dimension(x)))
{if(is.null(x$input))  NULL else  x$input[ , series, drop=FALSE]}

output.data.TSPADIdata <- function(x,series=seq(length=output.dimension(x)))
{if(is.null(x$output)) NULL else  x$output[ , series, drop=FALSE]}


#  default should work
# input.dimension.TSPADIdata <- function(x) {nseries( input.data(x))}
#output.dimension.TSPADIdata <- function(x) {nseries(output.data(x))}

# input.series.names, output.series.names default should work

identifiers.TSPADIdata <- function(obj) 
	{list(input=identifiers(obj$input), output=identifiers(obj$output))}
sourcedb.TSPADIdata <- function(obj) 
	{list(input=sourcedb(obj$input), output=sourcedb(obj$output))}
sourceserver.TSPADIdata <- function(obj) 
	{list(input=sourceserver(obj$input), output=sourceserver(obj$output))}
source.info.TSPADIdata <- function(obj) 
	{list(input=source.info(obj$input), output=source.info(obj$output))}


############################################################

#      Database interface for TSPADIdata  <<<<<<<<<<

############################################################



freeze.TSPADIdata <- function(x, timeout=60)
{ # This function retreives data from a PADI server using getpadi
  # See freeze.
  if (is.null(x$input))
    {z <- TSdata(output=freeze(x$output))
     z$source <- x
     return(z)
    }
  if (is.null(x$output)) 
    {z <-TSdata(input=freeze(x$input))
     z$source <- x
     return(z)
    }
  # now so that input and output are aligned ...

   if (! test.equal(attr(x$input,  "start") , attr(x$output, "start")))
         warning("input and output start values do no match. Using outputs.")
   if (! test.equal(attr(x$input,  "end") , attr(x$output, "end")))
         warning("input and output end values do no match. Using outputs.")
   if (! test.equal(attr(x$input,  "frequency"), attr(x$output, "frequency")))
        warning("input and output frequency values do no match. Using outputs.")

   if(attr(x$input, "pad.start") != attr(x$output, "pad.start") |
      attr(x$input, "pad.end")   != attr(x$output, "pad.end")   )
      warning ("input and output padding attibutes do not match. Using outputs")
      
   if(!(test.equal(attr(x$input, "use.tframe"),    attr(x$output, "use.tframe"))  &
        test.equal(attr(x$input, "start.server"),  attr(x$output, "start.server"))&
        test.equal(attr(x$input, "server.process"),attr(x$output, "server.process"))&
        test.equal(attr(x$input, "cleanup.script"),attr(x$output,"cleanup.script"))&
        test.equal(attr(x$input, "stop.on.error"), attr(x$output,"stop.on.error"))&
        test.equal(attr(x$input, "warn"),          attr(x$output, "warn") )))
      warning ("input and output server attibutes do not match. Using outputs")

  r <- freeze(modify(x$output, # output first so attributes are used
         append=list(series=x$input[1,],server=x$input[2,],db=x$input[3,],
	             transforms=x$input[4,],names=series.names(x$input))))
  r <- TSdata(input=r, output=r)
  r$source <- x
  input.data(r)  <-  input.data(r, series=ncol(x$output)+seq(length=ncol(x$input)))
  output.data(r) <- output.data(r, series=seq(length=ncol(x$output)))
  r
}


availability.TSPADIdata <- function(x, verbose=T, timeout=60)  
{# Indicate  dates for which data is available. 

 i <- if (0 ==  input.dimension(x)) NULL
      else availability( input.data(x), verbose=verbose)
 o <- if (0 == output.dimension(x)) NULL
      else availability(output.data(x), verbose=verbose)
 if (is.null(i) & is.null(o)) stop("No data.")
 invisible(list(start = rbind(i$start, o$start),
                end   = rbind(i$end, o$end),
                frequency=c(i$frequency, o$frequency),
                series=c(i$series, o$series)))
}


putpadi.TSdata   <- function (data, dbname, server=local.host.netname(), 
                   start.server=T, server.process=padi.server.process(), 
                   cleanup.script=padi.cleanup.script(),
                   series=series.names(data),
                   user=user.name(), passwd="",
                   stop.on.error=T, warn=T)   
  {#dbname and server can be a single string in which case it is applied to
   # all series. Otherwise it should be a structure like series: a list with
   # elements input and output, each vectors with a string for each series.

   # This function uses tfputpadi and returns an TSPADIdata object which can 
   #  be used to fetch the data from the database. tfputpadi in turn 
   #  uses putpadi.default.

   m <-input.dimension(data)
   p <-output.dimension(data)

   if(!is.list(dbname)) 
     {z <-dbname[1]
      dbname  <- list(input  = if (m==0) NULL else rep(z,m),
                      output = if (p==0) NULL else rep(z,p) )
     }

   if(!is.list(server)) 
     {z <-server[1]
      server <- list(input  = if (m==0) NULL else rep(z,m),
                     output = if (p==0) NULL else rep(z,p) )
     }

   if (m == 0) i <- NULL else
     {if(all (1 == start(input.data(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(input.data(data), list(start = start(input.data(data)), 
                frequency=frequency(input.data(data))))

      i <- tfputpadi(mat, server=server$input, dbname=dbname$input, 
         series=series$input,
         start.server = start.server, server.process = server.process, 
         cleanup.script = cleanup.script,
         user=user, passwd=passwd, stop.on.error=stop.on.error, warn=warn)   
     }
   if (p == 0) o <- NULL else
     {if(all (1 == start(output.data(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(output.data(data), list(start = start(output.data(data)),
                 frequency=frequency(data)))
      o <- tfputpadi(mat,  server=server$output, dbname=dbname$output, 
         series = series$output,
         start.server = start.server, server.process = server.process, 
         cleanup.script = cleanup.script,
         user=user, passwd=passwd, stop.on.error=stop.on.error, warn=warn)   

     }
   #This bypasses the constructor (structures are already built by tfputpadi):
   invisible(classed(list(input=i, output=o), c("TSPADIdata", "TSdata"))) # bypass constructor 
  }


#   The following function is supplied separately (with PADI ). The 
#   documentation is included here so it will integrate with DSE.



set.TSPADIdata <- function()
 {# prompt for input and output series identifiers, sets class, etc.
  cat("This function prompts for the names and database locations for\n")
  cat("input and output series, until an empty line is entered.\n")
  cat("If your model has no input or no output then return an empty line.\n\n")
  cat("Input (exogenous) variables...\n")
  i <- set.tfPADIdata(preamble=F)
  cat("Output (endogenous) variables...\n")
  o <- set.tfPADIdata(preamble=F)
  data <- classed(list(input=i, output=o),   # bypass constructor (set.TSPADIdata)
                  c("TSPADIdata", "TSdata"))  
  cat("The series may now be retrieved, in which case the data is\n")
  cat("  fixed as currently available, or they may be left `dynamic',\n")
  cat("  in which case they are retrieved using freeze.\n")
  cat("Retrieve data y/n:");key <- readline()
  if ((key =="y") | (key=="Y")) data <- freeze(data)
  data
}




retrieve.and.verify.data <- function(data.names,
             verification.data=verification.data, fuzz=1e-10)  
{# retrieve  data from a data base and do some verification. 
 # It is often useful if one of these is set:
 #    data.names$pad =T
 #    data.names$pad.end = T
 data <- freeze(data.names)
 #   check that data has not been rebased or otherwise messed up.
 if (0 != (input.dimension(data)))
   {s <-input.start(verification.data)
    e <-input.end(verification.data)
    error <- input.data(verification.data) -
               tfwindow(input.data(data),start=s, end=e, warn=F)
    if (fuzz < max(abs(error)) )
      {warning(paste("Retrieved input variables do not compare with the verification data.",
       "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
       " The maximum error is ", max(abs(error))))
       key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
       if (key=="y" | key=="Y")
         {z <- TSdata(input=error)
          tfplot(z, 
              select.inputs=(1:input.dimension(z))[apply(error,2,any)],
              select.outputs= 0)
         }
      key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
       if (key=="y" | key=="Y")
         {graph.data <- data
          output.data(graph.data) <-tfwindow(output.data(data),start=s,end=e)
          input.data(graph.data)  <-tfwindow(input.data(data), start=s,end=e)
          tfplot(verification.data, graph.data,
            select.inputs=(1:input.dimension(data))[apply(error,2,any)], select.outputs= 0)
         }
      }
   }
 s <-output.start(verification.data)
 e <-output.end(verification.data)
 error <-  output.data(verification.data) -
              tfwindow(output.data(data),start=s,end=e, warn=F)
 if (fuzz < max(abs(error))  )
   {warning(paste("Retrieved output variables do not compare with the verification data.",
    "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
    " The maximum error is ", max(abs(error))))
    key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
    if (key=="y" | key=="Y")
      {z <- TSdata(output=error)
       tfplot(z, select.inputs=0,
            select.outputs= (1:output.dimension(z))[apply(error,2,any)])
      }
    key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
    if (key=="y" | key=="Y")
      {graph.data <- data
       graph.output.data(data) <-tfwindow(output.data(data),start=s,end=e)
       graph.input.data(data)  <-tfwindow(input.data(data), start=s,end=e)
       tfplot(verification.data, graph.data, select.inputs=0,
            select.outputs= (1:output.dimension(data))[apply(error,2,any)])
      }
   }
data
}


#######################################################################

#    TS PADI interface tests (from Brief User's Guide)   <<<<<<<<<<

#######################################################################



TSPADI.function.tests <- function( verbose=T, synopsis=T,
      fuzz.small=1e-14, fuzz.large=1e-6, ets=F)
{# test for TSPADI access using simple.server
 # and if ets=T then run example from Brief User's guide (requires ets database)

 # These tests only check that the DSE structures work with PADI. For a more
 #   complete set of PADI tests see the file padi.s distributed 
 #   with the TS PADI software.


  if (synopsis & !verbose) cat("DSE TSPADI tests ...")

  scratch.db <-"zot123456.db"
  unlink(scratch.db)
  server <- local.host.netname()

 if (verbose) cat("DSE TSPADI test 0 ... ")
  if (check.padi.server(server))
     stop("A server is already running. Testing stopped. Use cleanup.padi.server() or kill.padi.server() to terminate it.")

  pid <- start.padi.server(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanup.padi.server(pid, cleanup.script="cleanup.simple.server"))

  # wait to ensure padi server is started
     for (i in 1:30)
       {if (check.padi.server(server)) break
        sleep(1)
       }

  exp1 <- tframed(matrix(1*exp(1:20),20,1), list(start=c(1950,1),freq=1))
#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  tframe(exp1) <- tframe(exp1)
  eg.put.data <- TSdata(input= exp1, 
                       output= tframed(tbind(2*exp1, 3*exp1),tframe(exp1)))
  series.names(eg.put.data) <- list(input="exp1", output=c("exp2","exp3"))

  if (any(input.series.names(eg.put.data) != "exp1"))
    stop("series.name setting is not working properly. Other tests will fail.")

  if (any(output.series.names(eg.put.data) != c("exp2","exp3")))
    stop("series.name setting is not working properly. Other tests will fail.")

#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  eg.put.data <- list(input= tsmatrix(exp1), 
#                      input.names="exp1",
#                      output= tsmatrix(2*exp1, 3*exp1), 
#                      output.names=c("exp2","exp3"))
  eg.names <- putpadi.TSdata(eg.put.data,
                      dbname=scratch.db, server=server,
                      start.server=T, server.process="simple.server", 
                      cleanup.script="cleanup.simple.server",
                      stop.on.error=T, warn=T )
  ok<-is.TSPADIdata(eg.names) 
  all.ok <- ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi server started\n")
    }

  if (verbose) cat("DSE TSPADI test 1 ... ")
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) & test.equal(eg.data, eg.put.data, fuzz=fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 2 ... ")

#If server= is supplied in the next, it should be "" and not NULL as previously 
eg.names <- TSPADIdata(input=c( "exp1","exp2"), output=c( "exp1","exp2","exp3"),
              frequency=1,
              db=scratch.db, stop.on.error=T, warn=T)

# z <- freeze(eg.names$input)
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) 
warning("skipping something broken")
#&
#    (max(abs(output.data(eg.data) - 
#              cbind(exp(1:20),2*exp(1:20),3*exp(1:20)) ))<fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 3 ... ")
  avail <- availability(eg.names, verbose=F)
  ok <- all(c(avail$start ==  t(matrix(c(1950,1),2,5)),
              avail$end   ==  t(matrix(c(1969,1),2,5)),
              avail$frequency ==  rep(1,5)))

  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  on.exit()
  cleanup.padi.server(pid, cleanup.script="cleanup.simple.server")

  if (synopsis) 
    {if (verbose) cat("All DSE TSPADI tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

if (ets)
{# test examples for TSPADI access (from Brief User's guide)
   # wait to ensure padi server is terminated
     for (i in 1:30)
       {if (!check.padi.server(server)) break
        sleep(1)
       }
   
  if (synopsis & !verbose) cat("DSE TSPADI/ets tests ...")

  if (verbose) cat("DSE TSPADI/ets test 1 ... ")
# this eventually does  getpadi("B1642", server="ets")

  eg2.DSE.data.names <- TSPADIdata(server="ets", db="",
#        output=c( "I37005"), output.names=c( "manuf.prod."), 
        output=c( "B1642"), output.names=c( "M1.sa."), 
        start.server=T, server.process="fame.server", 
        cleanup.script="cleanup.fame.server", stop.on.error=F, warn=T )

  z<- availability(eg2.DSE.data.names, verbose=F)
  if(any(z$end[,1]==1 & z$end[,2]==1))
    warning("Looks like some series have been discontinued.")

  eg2.DSE.data <- freeze(eg2.DSE.data.names)
  ok <- is.TSdata(eg2.DSE.data )
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI/ets test 2 ... ")
  eg3.DSE.data.names <- TSPADIdata(
     input ="B1642", input.transforms="percent.change", input.names="M1.sa.",
     output="B1650",output.transforms="percent.change", output.names="M2++.sa.",
     pad.start=F, pad.end =T,
     server="ets", db= "",
#        start.server=T, server.process="fame.server", 
#        cleanup.script="cleanup.fame.server",
     stop.on.error=F, warn=T )
  z<- availability(eg3.DSE.data.names, verbose=F)
  if(any(z$end[,1]==1 & z$end[,2]==1))
    warning("Looks like some series have been discontinued.")
 
  eg3.DSE.data <- freeze(eg3.DSE.data.names)
  ok <- is.TSdata(eg3.DSE.data )
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI/ets test 3 ... ")

  egJofF.1dec93.data.names <- TSPADIdata(
	input = c("B14017"), #etsmfacansim
	input.transforms= c("diff"),
	input.names=c("R90"),
#	output = c("P484549", "I37026", "b1627", "b14013",discont.
	output = c("B820600", "I37026", "b1627", "b14013",
		   "b4237", "D767608", "b3400", "M.BCPI", "M.JQIND", "M.CUSA0"),
                 # etscpi etsgdpfc etsmfacansim etsmfacansim etsdistscu
                 # etslabour etsforexch etsbcpi etsusa etsusa
  	output.transforms=c("percent.change", 
			"percent.change","percent.change",
			"diff", "diff", "percent.change",
			"percent.change", "percent.change",
			"percent.change", "percent.change"),
	output.names=c("CPI", "GDP", "M1", "RL", "TSE300", 
			"employment", "PFX", "com. price ind.", 
			"US ind. prod.", "US CPI"),
        server="ets", db= "",
#        start.server=T, server.process="fame.server", 
#        cleanup.script="cleanup.fame.server",
        stop.on.error=F, warn=T )

  z<- availability(egJofF.1dec93.data.names, verbose=F)
  if(any(z$end[,1]==1 & z$end[,2]==1))
    warning("Looks like some series have been discontinued.")
  egJofF.1dec93.data <- freeze(egJofF.1dec93.data.names)
  ok <- is.TSdata(egJofF.1dec93.data)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }


  if (synopsis) 
    {if (verbose) cat("All DSE TSPADI/ets tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }
}
invisible(all.ok)
}

#   2000/04/20 14:58:11 
###########################################################################

# Combination forecasting  functions.                       <<<<<<<<<<<<

###########################################################################


combine.and.forecast <- function(model, new.data,  
                      overlapping.period.forecast.tag="g", forecast.tag="f") 

{# model should be a TSmodel.
 # new data should be a list with $data and $overriding.data.
 # It can also contain elements data.tag and overriding.data.tag, character string
 #   tags which are passed along to construct.data.to.override.horizon.
 # $overriding.data is used in place of data and model forecasts to the horizon
 # for which it is available. $overriding.data should also include any input (policy)
 # variables to the forecast horzon.
 # best.guess in the result is a combination of available data, overriding.data,
 # and predictions. 
 # first splice and fill with model predictions.
 con.data <- construct.data.to.override.horizon(new.data, model, plot=F, 
                      forecast.tag=overlapping.period.forecast.tag) 
		      
 pred <-l(model, con.data ,predictT=dim(con.data$input)[1])$estimates$pred 
   # do residual analysis ?
# forecast<-forecast(l(model, con.data), percent=c(80,100,120), horizon=6, plot=F)
#   pchange<-percent.change(forecast[[1]],forecast[[2]],forecast[[3]], base=base,lag=12,cumulate=T,e=T)
 best.guess <-splice.tagged(con.data$output, pred, 
                  tag1=con.data$output.tags,tag2=forecast.tag) 
# the following result could also include con.data and pred, but it should be possible to
#    reconstruct these from the information in the list.

 invisible(list(model=model,
                data=new.data$data,
                overriding.data=new.data$overriding.data, 
                override=con.data$override,
                best.guess=best.guess))
}

reconstruct.combined.forecast <- function(combined.forecast) 
{# use the result of combine.and.forecast to re-do and verify results
 con.data <- construct.data.to.override.horizon(combined.forecast, combined.forecast$model, plot=F)
 pred <-l(combined.forecast$model, con.data ,predictT=dim(con.data$input)[1])$estimates$pred 
 best.guess <-splice.tagged(con.data$output, pred) 
 all(combined.forecast$best.guess==best.guess)
}

tfplot.combined.forecast <- function(combined.forecast,verbose=F, 
       start.=start(combined.forecast$data$output),
       Title="Projection", select.inputs=NULL, select.outputs=NULL, pause=T)
{# if verbose is T additional information is provided
 # if pause is true graphics are paused between pages.
   if (pause) dev.ask(T)
   if (verbose)
     {tfplot(combined.forecast$data, start.=start., Title="Data and combined.forecast")
      tfplot(combined.forecast$pred, start.=start.,
            Title="Model predictions (one step ahead for history)")
     }
   graph.data <- combined.forecast$data
   graph.data$output <- combined.forecast$best.guess
   if (is.null(select.inputs))  select.inputs  <- seq(dim(graph.data$input)[2])
   if (is.null(select.outputs)) select.outputs <- seq(dim(graph.data$output)[2])
   tfplot(graph.data, start.=start., Title="Projection", 
           select.inputs=select.inputs, select.outputs=select.outputs)
#   tfplot(combined.forecast$forecast[[2]],combined.forecast$forecast[[1]],
#         combined.forecast$forecast[[3]], start.=start.,
#         Title="Projection using future policy=most recent value and 20% higher and lower")
#   tfplot(combined.forecast$pchange[[2]],combined.forecast$pchange[[1]],
#         combined.forecast$pchange[[3]],start.=start., Title=
#    "Year over year percent change using future policy=most recent value and 20% higher and lower")
   invisible()
}

###########################################################################

# functions for misc. data retrieval, checking, and transformation <<<<<<<<<<<<

###########################################################################


construct.data.to.override.horizon <- function(new.data, model, plot=T, forecast.tag="f")
{# model should be a TSmodel.
 # new.data should be a list with $data and $overriding.data.
 # $overriding.data is used in place of $data and model forecasts to 
 # the horizon for which it is available. 
 #  Splice together $data and $overriding.data and if necessary
 #  calculate predictions for $overriding.data period and use them where $overriding.data
 #  or $data are not available, then return complete data set 
 #  to the end of the $overriding.data horizon, along with input data.
 #    (Typically the end of $overriding.data$output determines the periods
 #     for which forecast are combined and the end of $overriding.data$input
 #     determines how far into the future the model is used to extend the
 #     combined forecast. )
 #  Note that the $overriding.data is used in place of data in the 
 #  returned data set to allow for over-riding with anticipated data revisions.
 #  However, for any predictions during the combined.forecast period (ie. to augment
 #  $data and $overriding.data as returned by this function),  
 #  only $data is used and only to the last period for which observations
 #  for all variables are available.

 # if $overriding.data and $data overlap indicate override locations in 
 #     logical matrix dup:

 # tbind aligns the matrices
 dup <- tbind(output.data(new.data$data), output.data(new.data$overriding.data))
 if (!is.null(dup))
  {p <- output.dimension(new.data$data)
   dup <- (!is.na(dup[,1:p,drop=F])) & (!is.na(dup[,(p+1):(2*p),drop=F]))
  }

    # This can be used to provide a warning. eg
    #if (any(dup))
    #  {cat("WARNING:overriding is being used where data is available as follows:\n")
    #   print(dup)
    #  }

 z <- trim.na(output.data(new.data$data), end.=F)
 zz <- new.data$overriding.data$output
 z <- splice(zz,z)
 start. <- start(z)
 if (is.null(new.data$data$input)) z.in <-NULL
 else
   {# note that $overriding.data does not override actual data for $input, but 
    #  that could be changed by reversing the order in the next line. (extra  
    #  code is necessary to give a warning.)
    z.in <-trim.na(splice(input.data(new.data$data),
                          input.data(new.data$overriding.data)))
    start. <- latest.start(z, z.in)
    z.in <- tfwindow(z.in, start=start., warn=F)
    if (any(is.na(z.in)))
       stop(paste("Input (exogenous) series data cannot be specified as NA. (note ",
                  "differenced data requires an overlap of one period at the end of ",
                  "historical data and the beginning of monitoring overriding data.)"))
   }
 z <- tfwindow(z, start=start., warn=F)
 con.data <- TSdata(output=z,  input=z.in)

 # now augment $data and $overriding.data with model predictions for 
 #  the combined forecast period if necessary.
 if (any(is.na(output.data(con.data))))    
   {z <- TSdata(input = input.data(con.data),
                output= trim.na(output.data(new.data$data)))
    pred <- l(model,z, predictT= output.periods(con.data))$estimates$pred
    z <-splice.tagged(output.data(con.data),pred, 
                    tag1=con.data$output.tags, tag2=forecast.tag)
    output.data(con.data) <- z
   }

 con.data<- freeze(con.data)
 con.data$override <- dup
 if (plot && exists.graphics.device()) 
    {tfplot(con.data,start.=(end(output.data(data))-c(1,0)), 
           Title="Historical and overriding data data")
    }
  invisible(con.data)
}

get.overriding.data <- function(file="overriding.data", 
 first.input="",first.output="", second.output="", m=1, p=10)
{#Get other data eg(monitoring or other forecast data) 
  #   N.B. This cues on series names in the file
  # m is the number of input series
  # p is the number of output series
  z  <- dsescan(file=file,what=character())
  first.in   <- (1:length(z))[z==first.input] 
  if (0== length(first.in))
     stop(paste("Cannot find keying string:", first.input," in file", file))
  first.out  <- (1:length(z))[z==first.output] 
  if (0== length(first.out))
     stop(paste("Cannot find keying string:", first.output," in file", file))
  second.out <- (1:length(z))[z==second.output] 
  if (0== length(second.out))
     stop(paste("Cannot find keying string:", second.output," in file", file))
  input.periods <- (first.out-(first.in+m))/m     
  zz <- matrix(z[first.in:(first.out-1)],(input.periods+1),m)
  input.names <- zz[1,]
  input <-  matrix( as.numeric(zz[2:(1+input.periods),]), input.periods,m)
  dimnames(input) <- list(NULL,input.names)
  input <- tframed(input, list(start=as.integer(z[1:2]),frequency=12))
  output.periods<- second.out-(first.out+1)
  zz <- matrix(z[first.out:length(z)],(output.periods+1),p)
  output.names <- zz[1,]
  output <-  matrix( as.numeric(zz[2:(1+output.periods),]), output.periods,p)
  dimnames(output) <- list(NULL,output.names)
  output <- tframed(output, list(start=as.integer(z[1:2]),frequency=12))
  TSdata(input=input , output=output)
}


#tfplot.combined.forecast(combined.forecast,verbose=F, 
#      start.=start(combined.forecast$data$output),
#      Title="Projection", select.inputs=NULL, select.outputs=NULL)


restrict.overriding.data <- function(data, overriding.horizon=0)  
{#This function truncates overriding.data$output if it extends 
 #  overriding.horizon periods beyond the present. 
 year.mo <- c(date.parsed()$y,date.parsed()$m) - c(0,1) + c(0,overriding.horizon)
#check this - also note NAs should not be nec in overriding fame data
 data$output <-tfwindow(data$output, end=year.mo, warn=F )
 invisible(data)
}

###########################################################################

# functions for e-mail of results of combination forecasting <<<<<<<<<<<<

###########################################################################

combination.monitoring <- function(model, data.names,
   previous.data=NULL,
   overriding.data.names=NULL, 
   restrict.overriding.data=T, overriding.horizon=0,
   mail.list=NULL,
   error.mail.list=NULL,
   message.title="Combination Monitoring",
   message.subject="Combination Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables=series.names(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   future.input.data.tag="p",
   overriding.data.tag="m",
   overlapping.period.forecast.tag="g",
   forecast.tag="f",
   run.again=F,
   save.as=NULL)

{ # Step 0 - prepare message files and error checking
    error.message <- c(message.title, paste(date.parsed(), collapse="."),
              "An error condition occurred running combination.monitoring.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(mail(error.mail.list, subject=paste("error ", message.subject),
                 text= c(error.message, message)))
    if ( dseclass(model)[1] == "TSestModel" ) model <- model$model
    if (!is.null(data.names$pad.end))
       {if(!data.names$pad.end)
          warning("pad.end in data definition may disable retrieving all data.")
       } 
    else if (!is.null(data.names$pad))
       {if(!data.names$pad)
          warning("pad in data definition may disable retrieving all data.")
       } 

# The following line can be removed if the code works reliably
   mail(error.mail.list,subject=paste("checking ",message.subject),
                           text=paste(date.parsed(), collapse="."))

 # Step 1 - retrieve & check for updated data  or
 #            initialize system and if previous.data is NULL
    if (is.null(previous.data))
      {data <- freeze(data.names)
       message <- "Initializing combination monitoring:"   
       status <- "Combination monitoring initialized."   
      }
    else if (run.again)
      {data <-previous.data$data  
       overriding.update <- previous.data$overriding.data
       status <- "Combination monitoring re-run."   
      }
    else
      {updated.data<-check.for.value.changes(data.names,
                           verification.data=previous.data$data,
                           discard.current=T)
       if (is.null(overriding.data.names)) overriding.update<-list(F)
       else overriding.update<-check.for.value.changes(overriding.data.names,
                           verification.data=previous.data$overriding.data)
       if(updated.data[[1]] | overriding.update[[1]])
         {status <- "Combination monitoring updated."     
          if(updated.data[[1]])
            {data <-updated.data$data
             message <- c("data updates: ", 
                 series.names(data)$input[updated.data$input],
                 series.names(data)$output[updated.data$output])
            }
          if(overriding.update[[1]])
            {overriding.data <- overriding.update$data
             if(restrict.overriding.data & (!is.null(overriding.data$output))) 
                overriding.data <- restrict.overriding.data(overriding.data, 
                                 overriding.horizon=overriding.horizon)
             message <- c(message,"monitoring data updates: ",
             series.names(overriding.data)$input[ overriding.update$input],
             series.names(overriding.data)$output[overriding.update$output])
            }
         }
       else
         {on.exit()
          return(invisible(list(data=previous.data, 
                status="Combination monitoring updates not necessary.")))
         }
      }

 # Step 2 - check data and overriding data
   # Sometimes data is available as soon as there are any days in a month (with
   #   ignore on in Fame). The following 4 lines trim these, but that may not be
   #   the best way to handle them.
   year.mo <- c(date.parsed()$y,date.parsed()$m) - c(0,1)
   data  <- tfwindow(data,  end=year.mo, warn=F )
   fr <- c(frequency(data), 1)
      
   # check matching of starting date with end of available data.
   #   period for which all data is available in data
   end.ets <- end(trim.na(output.data(data))) 
   if (!is.null(overriding.data))
    {if (is.null(overriding.data$output))
     {overriding.data$output <- ts(matrix(NA, 1, output.dimension(data)),
                           end=end(data$output), 
                           frequency=frequency(data$output), 
                           names=dimnames(data$output)[[2]])
      if (!is.null(data$output.names))
         overriding.data$output.names <- data$output.names
     }
   else
     {if (!( (1+fr %*% end.ets) >= (fr %*%start(overriding.data$output))))
        stop(paste("Monitoring data (or NAs) must be indicated after ", end.ets))
      if (1== latest.end.index(output.data(data), output.data(overriding.data)))
         warning(paste("Overriding data file does not appear to be updated.",
         "True data is available past the end of the overriding data."))
    }}   

    if (is.null(overriding.data.names)) overriding.data <- NULL
    else
       overriding.data <- tagged(overriding.data,
          input.tags=future.input.data.tag, output.tags=overriding.data.tag)
    data <- tagged(data, input.tags=data.tag, output.tags=data.tag)

 # Step 3 - run forecast
   # warnings from this should be mailed!!!!
    combined.forecast<-combine.and.forecast(model, list(data, overriding.data),
           overlapping.period.forecast.tag=overlapping.period.forecast.tag, 
           forecast.tag=forecast.tag) 

 # Step 4 - write and mail files
    message <- c(message, "Projections are conditioned on forecast of ",
                            series.names(updated.data$data)$input, 
                          "                        with tranformation ",
                           data.names$input.transformations,
                          "The forecasts are now:")
    #starting and end period for plots & printing:
    start.<-(end(combined.forecast$data$output)+show.start) 
    end.  <-(end(combined.forecast$data$output)+show.end)
    # this can be too long if sufficient input data is not provided, so: 
    if ((fr %*% end(combined.forecast$best.guess)) < ((end.-c(0,1)) %*% fr))
       end.  <-end(combined.forecast$best.guess)

    report.variables$input<- 
            (report.variables$input == series.names(data.names)$input)
    report.variables$output<- 
            (report.variables$output == series.names(data.names)$output)


    rv <- tagged(
              combined.forecast$best.guess[,report.variables$output, drop=F],
              tags= (attr(combined.forecast$best.guess,"tags")
                             ) [,report.variables$output, drop=F])
    tframe(rv) <- tframe(combined.forecast$best.guess)
    inp <- splice(combined.forecast$data$input, 
                  combined.forecast$overriding.data$input,
                  tag1=data.tag, tag2=future.input.data.tag)
    rv <-tfwindow(cbind(inp,rv), start=start., end=end., warn=F) 
    message <- c(message,fprint(rv, digits=5, sub.title=data.sub.heading)) 

    if (any(combined.forecast$override))
       {message <- c(message, "WARNING: overriding data is being used where historical data is available as follows:",
              combined.forecast$override)
       }

#    print(tfwindow(tsmatrix(combined.forecast$data$input, combined.forecast$best.guess), 
#      start=start.), digits=print.digits)

# The following needs a postscipt viewer like gv or pageview
#    postscript(file=graphics.file, width=7, height=8, pointsize=14,
#        horizontal=F, onefile=F, print.it=F, append=F)
#    graph.combined.forecast(combined.forecast, start.=start.)
#    dev.off()
#    message <- c(message,"For graphic (in OpenWindows) type:\n    pageview ")
#    if("/" == substring(graphics.file,1,1) )
#             message <- c(message,graphics.file)
#    else
#      {pwd <- present.working.directory()
#       if("/tmp_mnt" == substring(pwd,1,8)) pwd <-substring(pwd,9)
#       message <- c(message,paste(pwd,"/",graphics.file, sep=""))
#      }
#    message <- c(message," in a command tool window. (Be patient. It takes a few seconds.)")

    if (!is.null(message.footnote)) message <-c(message, message.footnote)
    mail(mail.list, subject=message.subject, text= message)


 # Step 4 - clean-up
    if (!is.null(save.as)) 
      {assign(save.as, combined.forecast, where=1)
#       file.copy( graphics.file, save.as)   # save graph
      } 
    if (updated.data[[1]] ) previous.data$data   <- updated.data$data
    if ( overriding.update[[1]])
       previous.data$overriding.data<- overriding.update$data
    on.exit()
    #return latest data for comparison next time
    invisible(list(data=previous.data, status=status, message=message)) 
}


###########################################################################

# Tests function    <<<<<<<<<<<<

###########################################################################

combination.monitor.function.tests <- function( verbose=T, synopsis=T, 
         fuzz.small=1e-10,
         server.process = padi.server.process(),
         cleanup.script = padi.cleanup.script() )
{# Some of the tests here are really for functions defined in dse1 ... dse3
 #   but are not tested there to avoid assuming Fame access is available.
 # The main short coming of these tests is that they do not test
 #     functions which produce output or graphs.
 # These tests require access to Fame data bases and the files:
 #          monitoring.test.db    fake database 
 #          monitoring.test.info  comparison info. to check results
 #          monitoring.test.data  fake over-riding data 

 # Note also that the test data is not real data (it may have been differenced
 #  or otherwise transformed) and is only intended to test that functions
 #  work as originally specified. 

  server <- local.host.netname()
  db     <- paste(DSE.HOME,"/data/monitoring.test.db",sep="")

  if (synopsis & !verbose) cat("All combination monitor tests ...")
  all.ok <- T

  if (verbose) cat("combination monitor test 0 ... ")
  # simulated a database server
  pid <- start.padi.server(server=server, dbname=db, 
           server.process=server.process)
  on.exit(cleanup.padi.server(pid, cleanup.script=cleanup.script))

  # wait for server to start 
     for (i in 1:30)
       {if (check.padi.server(server)) break
        sleep(1)
       }
  ok <- T
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("combination monitor test 1 ... ")
  #  dbname=db would not be nec. with a public mode fame server

  test.data.names <- TSPADIdata(
      input="B14017", 
#      input.transforms= "diff",
       output=c( "P484549", "I37026", "lfsa201","b3400"), 
#      output.transforms= rep("percent.change",4),
      db=db, server=server,pad.end =T)

  source(paste(DSE.HOME,"/data/monitoring.test.info", sep=""))

  v.data <- verification.data
  v.data$output <- v.data$output[,c(1,2,6,7)]
  tframe(v.data$output) <- tframe(verification.data$output)
  ok <- is.TSdata(v.data)
  all.ok <- ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (loading verification data)\n")
    }

  if (verbose) cat("combination monitor test 2 ... ")
  data <-retrieve.and.verify.data(test.data.names, 
                                    verification.data=v.data)
  ok <- test.equal(data, ets.test.data, fuzz=fuzz.small)
  tags(data$input)  <- "data"
  tags(data$output) <- "data"
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (retrieve.and.verify.data)\n")
    }

  if (verbose) cat("combination monitor test 3 ... ")
  overriding.data <- get.overriding.data(
                   file=paste(DSE.HOME,"/data/monitoring.test.data", sep=""),
                   m=1, p=10,
                   first.input="diff(R90=B14017)", 
                   first.output="%change(CPI=P484549)",  
                   second.output="%change(GDP=I37026)"  )
  z.tf <-tframe(overriding.data$output)
  overriding.data$output <- overriding.data$output[,c(1,2,6,7)]
  tframe(overriding.data$output) <- z.tf
  tags(overriding.data$input) <- "over"
  tags(overriding.data$output) <- "over"
  ok <- test.equal(overriding.data, monitor.test.data, fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (get.overriding.data)\n")
    }

  if (verbose) cat("combination monitor test 4 ... ")
  combined.forecast<-combine.and.forecast(monitoring.test.model,
                  list(data=data, overriding.data=overriding.data)) 

#  ok <- fuzz.small > max(abs( combined.forecast$best.guess - 
#                            best.guess.test))
#  Rbug - kludge - the above chokes on - because classes are not the same and 
#  gives Error: invalid time series parameters specified

  ok <- fuzz.small > abs( sum(combined.forecast$best.guess) - 
                              sum(best.guess.test))
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (combine.and.forecast)\n")
    }

  if (synopsis) 
    {if (verbose) cat("All combination monitor tests completed")
     if (all.ok) cat(" OK\n\n")
     else    cat(", some FAILED!\n\n")
    }
  invisible(all.ok)
}


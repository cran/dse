#   2000/04/18 11:15:53  
###########################################################################

# Simple monitoring functions and data checking        <<<<<<<<<<<<

###########################################################################


check.for.value.changes <- function(data.names, verification.data,
     discard.current=F,
     ignore.before= NULL,
     fuzz=1e-10)
  { # Check if data is modified or if more data is available.
    # data.names is an object of class c("TSPADIdata","TSdata").
    # verification.data is an object of class TSdata.
    # T is returned for any series which has a modified value.
    #   NA in one series and not in the other is counted as modified.
    # If data are not the same length then series are padded with NA
    #  so new NAs on the end will not count as a change.
    # It is assumed that changes happen at the end (not the beginning) of
    #   the data. The data is trimmed at the beginning to the point where
    #   all series are available. (this simplifies padding the end with NA)
    # If ignore.before is not NULL it should indicate a year and period
    #   before which data is trimmed, no comparison is performed and the
    #   data before is not returned. If there are NAs at the beginning then
    #   trimming as described above may make the data even shorter than
    #   indicated by ignore.before.
    # discard.current controls whether current period data is considered.
    #  (some series are available for a period from the first day of the
    #   period, and are updated daily. Usually these should be discarded
    #   by setting discard.current=T)

   data <- freeze(data.names) 
   if (discard.current)
     {year.mo <- c(date.parsed()$y,date.parsed()$m) - c(0,1)
      data  <- tfwindow( data,  end=year.mo, warn=F )
     }
   if (!is.null(ignore.before)) 
     {data <- tfwindow(data, start= ignore.before)
      verification.data <-tfwindow(verification.data, start= ignore.before)
     }
   data <-trim.na(data, start.=T, end.=F)
   verification.data <-trim.na(verification.data, start.=T, end.=F)
   # which series are changed:
   if (is.null(input.series.names(data.names))) in.up <- NULL
   else
     {ld <-input.periods(data)
      lv <-input.periods(verification.data)
      l <- max(ld, lv)
      if (ld < l)
        input.data(data) <- ts(rbind(input.data(data),  
                                     matrix(NA,l-ld,input.dimension(data))),
                     start=start(input.data(data)),  frequency=frequency(data))
      if (lv < l)
        input.data(verification.data) <- ts(rbind(input.data(verification.data),
                                        matrix(NA,l-lv, input.dimension(data))),
                     start=start(input.data(verification.data)),
                     frequency=frequency(verification.data))
      z <- (is.na(input.data(data)) & is.na(input.data(verification.data)))   # both NA
    # next fixes an Splus bug (IMHO) that the dim is dropped for col matrix
      if (!is.matrix(z)) z <- array(z, dim(input.data(data)))
      z <- (abs(input.data(data) - input.data(verification.data)) <= fuzz) | z
      z <- z & !is.na(z)
      in.up <- !apply(z,2, all)
     }
   if (is.null(output.series.names(data.names))) out.up <- NULL
   else
     {ld <-output.periods(data)
      lv <-output.periods(verification.data)
      l <- max(ld, lv)
      if (ld < l)
        output.data(data) <- ts(rbind(output.data(data), 
                                      matrix(NA,l-ld, output.dimension(data))),
                         start=start(data), frequency=frequency(data))
      if (lv < l)
        output.data(verification.data) <- ts(
                                rbind(output.data(verification.data), 
                                      matrix(NA,l-lv, output.dimension(data))),
                     start=start(output.data(verification.data)),
                     frequency=frequency(verification.data))
      z <- ( is.na(output.data(data)) & is.na(output.data(verification.data)))    # both NA
    # next fixes an Splus bug (IMHO) that the dim is dropped for col matrix
      if (!is.matrix(z)) z <- array(z, dim(output.data(data)))
      z <- (abs(output.data(data) - output.data(verification.data)) <= fuzz) | z
      z <- z & !is.na(z)
      out.up <- !apply(z,2, all)
     }
   list(any(c(in.up,out.up)), input=in.up, output=out.up, data=data)   
  }

check.for.file.date.changes <- function(data.names, verification.dates)
  {# check file dates against previous dates
   # It is preferable to do file date checking with a Unix shell script rather 
   #   than in S, and then start S for further checks only when the time stamp
   #   on the database files has changed.
   up.in <-NULL
   if (!is.null(input.series.names(data.names)))
    {for (f in data.names$input$db) up.in <- c(up.in, file.date.info(f))
     inT <-any(verification.dates$input != up.in)
    }
   up.out <-NULL
   for (f in data.names$output$db) up.out <- c(up.out,file.date.info(f))
   outT <-any(verification.dates$output != up.out)
   list( any(c(inT,outT)), input=inT, output=outT, 
         dates=list(input=up.in, output=up.out))
  }



simple.monitoring <- function(model, data.names, 
   previous.data=NULL,
   mail.list=NULL,
   error.mail.list=user.name(),
   message.title="Simple Monitoring",
   message.subject="Simple Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables= series.names(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   forecast.tag="f",
   run.again=F,
   save.as=NULL)

{# Step 0 -  prepare message files and error checking
    error.message <- c(message.title, paste(date.parsed(), collapse="."),
              "An error condition occurred running simple.monitoring.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(mail(error.mail.list,
                 subject=paste("error ",message.subject),
                 text= c(error.message, message)))
    if ( dseclass(model)[1] == "TSestModel" ) model <- TSmodel(model)
    if (!is.null(data.names$pad.end))
       {if(!data.names$pad.end)
          warning("pad.end in data definition may disable retrieving all data.")
       } 
    else if (!is.null(data.names$pad))
       {if(!data.names$pad)
          warning("pad in data definition may disable retrieving all data.")
       } 

# The following line is useful for debugging
#mail(error.mail.list, subject=paste("checking ",message.subject), 
#                         text=paste(date.parsed(), collapse="."))

 # Step 1 - retrieve & check for updated data  or
 #            initialize system and if previous.data is NULL
    if (is.null(previous.data))
      {data <- freeze(data.names)
       message <- "Initializing simple monitoring:"   
       status <- "Simple monitoring initialized."   
      }
    else if (run.again)
      {data <-previous.data  
       status <- "Simple monitoring re-run."   
      }
    else
      {updated.data<-check.for.value.changes(data.names,
                           verification.data=previous.data,
                           discard.current=T)
       if(updated.data[[1]])
         {data <-updated.data$data
          message <- c("data updates: ", 
               input.series.names(data)[ input.data(updated.data)],
              output.series.names(data)[output.data(updated.data)])
          status <- "Simple monitoring updated."   
         }
       else
         {on.exit()
          return(invisible(list(data=previous.data, 
                status="Simple monitoring updates not necessary.")))
         }
      }

 # Step 2 - check data
   # Sometimes data is available as soon as there are any days in a month (with
   #   ignore on in Fame). The following 4 lines trim these, but that may not be
   #   the best way to handle them.
   year.mo <- c(date.parsed()$y,date.parsed()$m) - c(0,1)
   data  <- tfwindow(data,  end=year.mo, warn=F )

 # Step 3 - run forecast
   pred<-forecast(model, data)$forecast[[1]]
   pred <-splice.tagged(output.data(data), pred, tag1=data.tag,tag2=forecast.tag) 
 
 # Step 4 - generate report and mail
    message <-c(message,"The forecasts are now:")
    #starting and end period for plots & printing:
    start.<-(output.end(data)+show.start) 
    end.  <-(output.end(data)+show.end)

    report.variables$input<- 
            (report.variables$input == input.series.names(data.names))
    report.variables$output<- 
            (report.variables$output == output.series.names(data.names))
    rv <- tagged(pred[,report.variables$output, drop=F],
                 tags= (attr(pred,"tags")) [,report.variables$output, drop=F])
    tframe(rv) <- tframe(pred)
    inp <-tagged(input.data(data)[,report.variables$input, drop=F],tags= data.tag)

    tframe(inp) <-  tframe(input.data(data))

    rv <- tfwindow( tbind( inp, rv), start=start., end=end., warn=F)   
    message <- c(message,fprint(rv, digits=5, sub.title=data.sub.heading)) 

    if (!is.null(message.footnote)) message <-c(message, message.footnote)
    mail(mail.list, subject=message.subject, text= message)

 # Step 4 - clean-up
    if (!is.null(save.as)) 
       assign(save.as,list(model=model, data=data, pred=pred), where=1)
    on.exit()
    #return latest data for comparison next time. Note - the forecast is NOT
    # returned (but may be saved above).
    invisible(list(data=data, status=status, message=message)) 
}

watch.data <- function(data.names, 
   previous.data=NULL,
   mail.list="gilp",
   error.mail.list=NULL,
   message.title="Data Monitor\n",
   message.subject="Data Monitor",
   message.footnote=NULL)

{# monitors data bases and check series for changes with e-mail of results.
 # this should be used with a script which watches for file date changes.
 #  ( see example in the file watch.data.readme)
 # data.names is a TSdata (names) object.
 # mail.list and error.mail.list should be single strings (not vectors)
 # If mail.list is null then mail is not sent (useful for testing).
 #  but the string can contain multiple user ids for mail
 # previous.data must normally be supplied. If it is not (ie. if it is NULL)
 # then the system will be initiallized and the returned result will be
 # the previous.data for the next time the function is called.

 # Step 0 - prepare message files 
    error.message <- c(message.title, paste(date.parsed(), collapse="."),
              "An error condition occurred running watch.data.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(mail(error.mail.list, subject=paste("error ",message.subject),
                 text= c(error.message, message)))

 # Step 1 - retrieve & check for updated data 

    data.names <- modify.TSPADIdata(data.names, pad.end=T)
    #  Initialize system and exit if previous.data is NULL
    if (is.null(previous.data))
      {current.data <- freeze(data.names)
       on.exit()
       #return latest data for comparison next time
       return(invisible(list(data=current.data,
           status="System watch.data initialized."))) 
      }
    update<-check.for.value.changes(data.names,
                           verification.data=previous.data$data,
                           discard.current=F)
    if (!update[[1]] )
        {on.exit()
         return(invisible(list(data=previous.data$data, 
             status="No data updates.")))
        }
    else
       message <- c(message, "data updates: ", 
              output.series.names(update$data)[update$output],)

 # Step 2 - mail 
    if(!is.null(message.footnote)) message <- c(message,message.footnote)
    mail(mail.list, subject=message.subject, text= message)

 # Step 3 - clean-up
    on.exit()
    #return latest data for comparison next time
    invisible(list(data=update$data, status="Data has been updated.")) 
}


###########################################################################

# Tests function for data retrieval for simple monitoring    <<<<<<<<<<<<

###########################################################################


simple.monitor.function.tests <- function( verbose=T, synopsis=T, 
         fuzz.small=1e-14, fuzz.large=1e-8,
         server.process = padi.server.process(),
         cleanup.script = padi.cleanup.script() )
{# Some of the tests here are really for functions defined in dse1 ... dse3
 #   but are not tested there to avoid assuming TSPADI (or Fame) access is
 # available. The main short coming of these tests is that they do not test
 #     functions which produce output or graphs.
 # These tests require access to Fame data bases and the files:
 #          monitoring.test.db    fake database 
 #          monitoring.test.info  comparison info. to check results

 # Note also that the test data is not real data (it may have been differenced
 #  or otherwise transformed) and is only intended to test that functions
 #  work as originally specified. 

  server <- local.host.netname()
  db     <- paste(DSE.HOME,"/data/monitoring.test.db",sep="")

  if (synopsis & !verbose) cat("All simple monitor tests ...")
  all.ok <- T

  if (verbose) cat("simple monitor test 0 ... ")
  # simulate a database server
  pid <- start.padi.server(server=server,
           dbname=db, 
           server.process=server.process)
  on.exit(cleanup.padi.server(pid, cleanup.script=cleanup.script))

  # wait for server to start 
     for (i in 1:30)
       {if (check.padi.server(server)) break
        sleep(1)
       }
  ok <- T
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 1 ... ")
  #  db=db would not be nec. with a public mode fame server   
  test.data.names <- TSPADIdata(
      input  ="B14017", 
      output = c( "P484549", "I37026", "lfsa201","b3400"), 
      server=server, db=db, pad.end =T)
   
  z <-availability(test.data.names, verbose=F) 
  ok <- all(c(z$start == t(matrix(c(1974,2),2,5)), 
              z$end   == t(matrix(c(1993,9),2,5)), 
              z$freq==rep(12,5) ))
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


# the following sets ets.test.data, monitor.test.data, verification.data
#      and  monitoring.test
  source(paste(DSE.HOME,"/data/monitoring.test.info", sep=""))

  if (verbose) cat("simple monitor test 2 ... ") 
  v.data <- verification.data
  output.data(v.data) <- output.data(v.data)[,c(1,2,6,7)]
  tframe(output.data(v.data)) <- tframe(output.data(verification.data))
  ok <- is.TSdata(v.data)
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("simple monitor test 3 ... ")
  hist.data <-retrieve.and.verify.data(test.data.names, 
                                    verification.data=v.data)
  ok <- test.equal(hist.data, ets.test.data, fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 4 ... ")
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
                   previous.data=NULL, mail.list=NULL, error.mail.list=NULL) 
  ok <-  monitoring$status == "Simple monitoring initialized."   
  if (verbose) cat("\n This test produces a warning: Input is not longer than output data. No forecasts produced...")
  # note that the following does not result in forecasts (and the forecast
  #   function produces a warning) because the input data does not extend
  #   beyond the output data.
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
           previous.data=monitoring$data, mail.list=NULL, error.mail.list=NULL) 
  ok <- ok & (monitoring$status == "Simple monitoring updates not necessary.")
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
                 previous.data=monitoring$data, 
                 mail.list=NULL, error.mail.list=NULL, run.again=T) 
  ok <- ok & (monitoring$status == "Simple monitoring re-run.")
  ok <- ok & monitoring$message[7] == 
          "1993 Sep   0.110000   0.383440   0.397520   0.355500   0.947460 "
  ok <- ok & sum(output.data(monitoring$data))==235.64806565791809589
  output.data(monitoring$data) <- 
               tfwindow(output.data(monitoring$data), end=c(1993,8))
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
          previous.data=monitoring$data, mail.list=NULL, error.mail.list=NULL) 
  ok <- ok & (monitoring$status == "Simple monitoring updated.") &
      sum(output.data(monitoring$data)) == 235.64806565791809589
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 5 ... ")

  watch <- watch.data(test.data.names, previous.data=NULL, mail.list=NULL)
  ok <- (watch$status == "System watch.data initialized.") & 
         sum(output.data(watch$data))== 235.64806565791809589
  watch <- watch.data(test.data.names, previous.data=watch, mail.list=NULL)
  ok <- ok & (watch$status == "No data updates.") & 
           sum(input.data(watch$data))== -4.1300000572204575988
  watch$data <- tfwindow(watch$data, end=c(1993, 8))
  watch <- watch.data(test.data.names, previous.data=watch, mail.list=NULL)
  ok <- ok & (watch$status == "Data has been updated.") & 
          sum(output.data(watch$data))== 235.64806565791809589

  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All simple monitor tests completed")
     if (all.ok) cat(" OK\n\n") else    cat(", some FAILED!\n\n")
    }
  invisible(all.ok)
}

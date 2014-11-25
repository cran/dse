
#  Functions in this file are defined for S (and differently in R).

is.S <- function(){is.Svanilla() | is.Splus() }
is.Svanilla <- function(){!exists("version")}
is.Splus <- function(){exists("version") && is.null(version$language)}
is.Splus.pre3.3 <- function()
   { ## <= 3.2
    if( !exists("version")) 
         { #-- `Vanilla' S (i.e. here "S version 4")
           #- this now works for  S version 4  (this is not S-plus 4.0, maybe 
           #             part of S-plus 5.0 !):
           lv <- nchar(Sv <- Sversion())
           r <- list(
	      major = substring(Sv, 1,1),
	      minor = substring(Sv, lv,lv))
         }
   else 
     {r <- version
      r$minor <- as.numeric(r$minor)
      r$major <- as.numeric(r$major)
     }
    is.Splus() &&  ((r$major+.1*r$minor) < 3.3)
   }


   is.R <- function() {F}

    # faking R compatibility:
    
    require <- function(...){T} 
    Sys.getenv <- getenv

    Platform <- function()
      {list( OS.type = if(charmatch("MS Windows", version$os, nomatch=0)) "windows" else
       if(exists("unix")) "unix"  else "unkown")
      }
   
    Sys.info <- function()
     {if( !exists("version")) 
         { #-- `Vanilla' S (i.e. here "S version 4")
           #- this now works for  S version 4  (this is not S-plus 4.0, maybe 
           #             part of S-plus 5.0 !):
           lv <- nchar(Sv <- Sversion())
           r <- list(
	      major = substring(Sv, 1,1),
	      minor = substring(Sv, lv,lv))
         }
       else 
         {r <- version
          r$minor <- as.numeric(r$minor)
          r$major <- as.numeric(r$major)
         }
       if      (is.Splus())    r$language <- "default"
       else if (is.Svanilla()) r$language <- "S"
       r$OSversion <- paste(system("uname -s", intern=T), 
           system("uname -r | sed -e 's/\\.\.\*//'", intern=T), sep="")
       r$OStype    <- Platform()$OS.type
       r$nodename  <- system("uname -n", intern=T)
       r$login     <- system("whoami", intern=T) # not always correct
       if(0 == length(r$login)) r$login <- "unknown"
       r$user      <- Sys.getenv("USER")
       if(0 == length(r$user)) r$user <- "unknown"
       class(r) <- "Sys.info"
       r
      }
      

    .Syskern.options <- list(mailer = if ( "unix" == Platform()$OS.type )
                               "mailx" else "none" )
    getOption <- function(x) {.Syskern.options[[x]]}

   
    .SPAWN <- TRUE

    syskern.rm <- unlink

     

  # in S these only work in Unix so far.
    if("unix" == Platform()$OS.type) {
      system <- function(cmd, intern=F) {unix(cmd, output.to.S = intern)}
                 
      file.exists <- function(file) {1 == unix(paste("if [ -f ", file, 
           " ] ; then (echo 1) ; else (echo 0); fi"))}
      date.parsed <- function() 
        {d <- parse(text=system("date '+%Y %m %d %H %M %S'", intern=T),white=T)
         list(y=  eval(d[1]),
              m=eval(d[2]),
              d= eval(d[3]),
              H= eval(d[4]),
              M= eval(d[5]),
              S= eval(d[6]),
              tz=system("date '+%Z'", intern=T))
        }

      getwd <- function(){system("pwd", intern=T)} 
      Sys.sleep <- function(n) {system(paste("sleep ", n))} 
      file.copy <- function(from, to)system(paste("cp ", from, to)) 
      file.date.info <- function(file.name)
     	{# This could be a lot better. It will fail for files older than a year.
      	# Also, a returned format like date() below would be better.
      	mo <- (1:12)[c("Jan","Feb","Mar","Apr","May", "Jun","Jul","Aug", "Sep",
         "Oct","Nov","Dec") ==
          substring(system(paste("ls -l ",file), intern=T),33,35)]
      	day <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),37,38))
      	hr  <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),40,41))
      	sec <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),43,44))
      	c(mo,day,hr,sec)
     	}
      }

    if("windows" == Platform()$OS.type)
      {file.date.info <- Sys.sleep <- system <-
         getwd <- file.copy <- function(cmd) 
         {stop("syskern system calls not implemented for this operating system.")}
      }
    

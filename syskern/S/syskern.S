
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
    
    require <- function(package, ...)
      {if(package == "padi")
           {if("" == Sys.getenv("PADI_LDLIB")) return(FALSE)
	    else if (! exists("load.padi"))    return(FALSE)
	    load.padi()  
	    return(TRUE)}
       else TRUE 
      }
        
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
    

if (is.S())
  {
    RNGkind <- function(kind=NULL, normal.kind=NULL)
      {# With a null argument this returns the current kind and normal.kind.
       # kind = "default" resets the RNG to the Splus default.
       # Splus does not allow arbitrary setting of .Random.seed so
       #     .RandomSeed is used. The existence of .RandomSeed is
       #     used to indicate whether an alternate RNG is used and the
       #     first element of .RandomSeed indicates the generator.
       #  Note this does not always work with For loops if where=0 is used.
       old <- if(exists(".RandomSeed"))  
       		c("Wichmann-Hill")[1+.RandomSeed[1]] else "default"
       if ( !is.null(kind))
         {# set the RNG kind
          if (kind == "default") 
            {if(exists(".RandomSeed", where=1)) remove(".RandomSeed", where=1)}
          else if (kind == "Wichmann-Hill") 
             assign(".RandomSeed", c(0, as.integer(100000*runif(3))), where=1)
          else stop("Only Wichmann-Hill, default or NULL supported for kind.")
         }
       old.normal <- if(exists(".RNORMtransform", where=1))
       		 .RNORMtransform else  "default"
       if ( !is.null(normal.kind)) 
          {if(exists(".BM.seed", where=1)) remove(".BM.seed", where=1)
	   if (normal.kind == "Box-Muller")
	         assign(".RNORMtransform", normal.kind, where=1)
           else if (normal.kind == "default")  
              {if(exists(".RNORMtransform", where=1))
                  remove(".RNORMtransform", where=1)
	       if(exists(".RandomSeed", where=1))
	         {warning("kind also set to default as required by default normal.kind") 
		  remove(".RandomSeed", where=1)
		 } 
	      }
           else stop("Only Box-Muller, default or NULL supported for normal.kind.")
          }
       c(old, old.normal)
      }

    if (!exists("set.seed.Splus")) set.seed.Splus <- set.seed

    set.seed <- function(seed=NULL)
      {# with a null argument this also serves as get.seed.
       kind <- RNGkind()
       if ( is.null(seed)) 
         {if (kind[1] == "default") seed <-.Random.seed
          else                      seed <-.RandomSeed[-1]
         }
       else
         {# set seed
          if (kind[1] == "default") 
             {if (1==length(seed)) set.seed.Splus(seed)
              else                 assign(".Random.seed", seed, where=1)#default
             }
          else if (kind[1] == "Wichmann-Hill") 
             {if (3 != length(seed))
                 stop("seed length is not consistent with kind Wichmann-Hill.")
              #Note this does not always work with For loops if where=0 is used.
              assign(".RandomSeed", c(0,seed), where=1)
             }
          else stop("seed does not match RNG kind.")
         }
       seed
      }
 
 
    set.RNG <- function(kind=NULL, seed=NULL, normal.kind=NULL)
      {# with a null argument this also serves as get.RNG 
        old <- list(kind=RNGkind()[1], normal.kind=RNGkind()[2],
	             seed=set.seed())
        if (is.null(kind) & is.null(seed) & is.null(normal.kind)) return (old)
	if (is.list(kind)) 
          {seed        <- kind$seed
	   normal.kind <- kind$normal.kind
	   kind        <- kind$kind
	  }
	RNGkind(kind=kind, normal.kind=normal.kind)
	set.seed(seed)
	old
      }


    if (!exists("runif.default")) runif.default <- runif
    runif <- function(n, min=0, max=1)
       {# This typically just calls runif.default, but allows using other
        # RNGs to generate the same sequence in R and S.
        # eg: set.RNG(seed=c(1:3), kind="Wichmann-Hill")
        #     runif(10)

        if(RNGkind()[1] == "default")  return(runif.default(n, min=min, max=max))
        else seed <- set.seed() # returns the setting
        kind <-  RNGkind()[1]
        if(kind == "Wichmann-Hill")
           {out <- numeric(n)  
            if (3 != length(seed)) stop("seed setting is not consistent with RNG.")
            x <- seed[1]; y <- seed[2]; z <- seed[3]
            for(i in 1:length(out))
               {x <- (171*x) %% 30269
                y <- (172*y) %% 30307
                z <- (170*z) %% 30323
                out[i] <- (x/30269 + y/30307 + z/30323) %% 1.0
               }
            set.seed( c(x,y,z))
           }
        else stop("runif RNG kind not supported.")
        out
       }


if (!exists("rnorm.default")) rnorm.default <- rnorm

rnorm <- function(n, mean=0, sd=1, compiled=F)
   {# This typically just calls rnorm.default, but provides the possibility of 
    # using Wichmann-Hill to generate the same runif sequence in R and S and 
    #    then generate the same normally distributed numbers with Box-Muller.
    # eg: set.RNG(seed=1:3, kind="Wichmann-Hill", normal.kind="Box-Muller")
    #   where 1:3 should be a valid seed.
    # This replicates R values, given by
    #   set.RNG(seed=1:3, kind="Wichmann-Hill", normal.kind="Box-Muller"),


    if(RNGkind()[2] != "Box-Muller") return(rnorm.default(n, mean=mean, sd=sd))
    else
      {if(n==0) return(numeric(0))
#       if(exists(".BM.seed", envir=.GlobalEnv)) 
       if(exists(".BM.seed", where=0)) 
         {out <- get(".BM.seed", where=0)
	  remove(".BM.seed", where=0)
	 }
       else out <- NULL
       # next should be true except when n==1 and an odd value has been saved
       if (length(out) < n) 
         {rv <- runif(n-length(out) + (n-length(out))%%2)
          rv <- matrix(rv, 2, length(rv)/2)
          rv <- c( rbind(sqrt(-2*log(rv[2,])) * cos(2*pi*rv[1,]),
                         sqrt(-2*log(rv[2,])) * sin(2*pi*rv[1,])))
          out <- c(out, rv)
	 }
       if (1 == (length(out) - n)) 
          {#drop last point and keep for next call
	   assign(".BM.seed", out[length(out)], where=0)
	   out <- out[-length(out)]
	  }
       if(n !=length(out)) stop("something is rotten in the state of rnorm.")
      }
    mean + out*sd
   }

  set.RNG(kind="default", normal.kind="default")

  }   # end of if is.S

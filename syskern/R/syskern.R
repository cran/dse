
system.info <- function()
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
   r$OSversion <- OSversion()
   r$OStype    <- OStype()
   r
  }


###########################################################

#    1/  Functions are for identifying S or R and flavours.

###########################################################

#Note:  It is tempting to use system.info as defined above to define the 
#       following, but there is a bootstrapping problem to 
#       solve (system.info uses is.XXX() so these cannot rely on it).

if (! exists("is.R"))
 {is.R <- function()
     {exists("version") && !is.null(vl <- version$language) && vl == "R" }
 }

is.R.pre0.60 <- function()
  {is.R() && ((as.numeric(version$major)+.01*as.numeric(version$minor)) <0.60) }
is.R.pre0.63.2 <- function()
  {is.R() && ((as.numeric(version$major)+.01*as.numeric(version$minor)) <0.623)}


is.S <- function(){is.Svanilla() | is.Splus() }
is.Svanilla <- function(){!exists("version")}
is.Splus <- function(){exists("version") && is.null(version$language)}
is.Splus.pre3.3 <- function()
   { ## <= 3.2
    is.Splus() &&  ((system.info()$major+.1*system.info()$minor) < 3.3)
   }

###########################################################

#    2/  Functions are for identifying the operating system.

###########################################################

if (is.R())
   {OStype <- function()
      {if ("windows" == .Platform$OS.type) return("MS Windows")
       else if("mac" == .Platform$OS.type) return("Macintosh") #needs to be checked
       else if("unix"== .Platform$OS.type) return("Unix")
      }
   }

if (is.S())
   {OStype <- function()
      {if(charmatch("MS Windows", version$os, nomatch=0))
                                return("MS Windows")
       else if(charmatch("Macintosh",  version$os, nomatch=0))
                                return("Macintosh") # needs to be checked
       else if(exists("unix"))  return ("Unix") 
      }
   }


is.MSwindows <- function(){OStype() == "MS Windows"}
is.Mac       <- function(){OStype() == "Macintosh" }  
is.unix      <- function(){OStype() == "Unix" }  

{
if (is.unix())
  {OSversion <- function()
    {paste(system.call("uname -s"), 
           system.call("uname -r | sed -e 's/\\.\.\*//'"), sep="") }
  }
else if(is.MSwindows())
  {if (is.R())
     {OSversion <- function() 
        {# This is not great since NT is not distinguished but
         #    is.Win32() below will work ok
         if("Win32"== machine()) return("MS Windows 95")
         else return ("unkown")
        }
     }
   if (is.S())
     {OSversion <- function() 
        {if("MS Windows 3.1"==version$os) return("MS Windows 3.1")
         if("MS Windows 95" ==version$os) return("MS Windows 95")
         if("MS Windows 98" ==version$os) return("MS Windows 98")
         if("MS Windows NT" ==version$os) return("MS Windows NT")
         else return ("unkown")
        }
     }
  }
else OSversion <- function() "unknown"  
}


# Other is.XXX() should be added here.

# determining Unix flavours doesn't seem to be too important but ...
is.Sun4 <- function() {is.unix() && OSversion() == "SunOS4" }
is.Sun5 <- function() {is.unix() && OSversion() == "SunOS5" }
is.Linux <- function(){is.unix() && OSversion() == "linux"} 

# Windows flavours may be more important but these are untested !!!
is.Win3.1 <- function(){is.MSwindows() && OSversion() == "MS Windows 3.1"} 
is.Win95  <- function(){is.MSwindows() && OSversion() == "MS Windows 95"} 
is.WinNT  <- function(){is.MSwindows() && OSversion() == "MS Windows NT"} 
is.Win32  <- function(){is.Win95() | is.WinNT() } 





###########################################################

#    3/  Functions depending only on the 
#         differences between S and R

###########################################################

if(is.S())
   {if(is.unix())system.call  <- unix   
    global.assign <- function(name, value) {assign(name,value, where = 1)}
    .SPAWN <- TRUE
    exists.graphics.device <- function(){dev.cur() !=1 }
    open.graphics.device  <- function(display=getenv("DISPLAY"))
                                 {openlook(display) }
    #                            {motif(display) }


    "list.add<-" <- function(x, replace, value)
       {# replace or add elements to a list.
        x[replace] <- value
        # x[[replace]] <- value  would be more logical but doesn't work
        x
       }
    require <- function(...){T} # faking R compatibility
   }
        
if(is.R()) 
   {#tempfile <- function(f)
    #   {# Requires C code also from Friedrich Leisch not in version 0.15 of R.
    #    d<-"This is simply a string long enough to hold the name of a tmpfile";
    #     .C("tmpf", as.character(d))[[1]]
    #    }

    if (is.R.pre0.60())
        {tempfile <- function(pattern = "file") 
                {system(paste("for p in", paste(pattern, collapse = " "), ";",
                       "do echo /tmp/$p$$; done"),
                 intern = TRUE)
                }
        }

#  no longer needed  unlink <- function(file) system.call(paste("rm -fr ", file))
    global.assign <- function(name, value) 
                          {assign(name,value, envir=.GlobalEnv)}
    synchronize <- function(x){NULL} # perhaps this should do something?
    .SPAWN <- FALSE
    dev.ask <- function(ask=T){par(ask=ask)}
    if (is.R.pre0.63.2())
         exists.graphics.device <- function(){exists(".Device")}  
    else exists.graphics.device <- function(){dev.cur() !=1 }
    open.graphics.device  <- function(display=getenv("DISPLAY")) {x11(display)}

   "list.add<-" <- function(x, replace, value)
     {# replace or add elements to a list. 
      if (is.numeric(replace))
        {# x<- do.call("default.[[<-", list(x,replace,value))   # use default
         x[[replace]] <- value
         return(x)
        }
      if (is.null(value))  value <- list(NULL)
      if (!is.list(value)) value <- list(value)
      if (1 == length(value)) 
       {for (i in seq(length(replace)))
          x<- do.call("$<-", list(x,replace[i],value[[1]]))
       }
      else
        {if(length(value) != length(replace) )
         stop("number of replacement values != number of elements to replace")
         for (i in seq(length(replace)))
            x<- do.call("$<-", list(x,replace[i],value[[i]]))
        }
      x
     }
 }


###########################################################

#    4/  Functions depending only on the 
#         differences among operating system.

###########################################################

if(is.unix())
  {sleep <- function(n) {system.call(paste("sleep ", n))} # pause for n seconds
   present.working.directory <- function(){system.call("pwd")} #present directory
   whoami <- function(){system.call("whoami")} # return user id (for mail)
   local.host.netname <- function() {system.call("uname -n")}

   mail <- function(to, subject="", text="")
     {# If to is null then mail is not sent (useful for testing).
      file <- tempfile()
      write(text, file=file)
      if(!is.null(to))
         system.call(paste("cat ",file, " | mail  -s '", subject, "' ", to))
      unlink(file)
      invisible()
     }

   file.copy <- function(from, to)system.call(paste("cp ", from, to)) #copy file

   file.date.info <- function(file.name)
     {# This could be a lot better. It will fail for files older than a year.
      # Also, a returned format like date() below would be better.
      mo <- (1:12)[c("Jan","Feb","Mar","Apr","May", "Jun","Jul","Aug", "Sep",
         "Oct","Nov","Dec") ==
          substring(system.call(paste("ls -l ",file)),33,35)]
      day <- as.integer(substring(system.call(paste("ls -l ",file.name)),37,38))
      hr  <- as.integer(substring(system.call(paste("ls -l ",file.name)),40,41))
      sec <- as.integer(substring(system.call(paste("ls -l ",file.name)),43,44))
      c(mo,day,hr,sec)
     }

}

if(is.MSwindows())
  {system.call  <- function(cmd) 
         {stop("system calls must be modified for this operating system.")}
   sleep <- system.call 
   present.working.directory <- system.call
   whoami <- system.call
   file.copy <- system.call
   file.date.info <- system.call
  }



###########################################################

#    5/  Functions depending on both R/S and the 
#         differences among operating system.

###########################################################

if(is.unix())
  {if(is.R()) 
     {#unix <- function(cmd) system(cmd, intern=T)
      # unix() is now a function in R but deprecated in favour of system()
      # (This is a bit dangerous, as these calls may be system dependent.)

      system.call <- function(cmd) system(cmd, intern=T)

  # the following date function might be made system independent as a C call.
      date.parsed <- function() 
        {d<-parse(text=strsplit(
              system.call("date \'+%Y %m %d %H %M %S\'")," ")[[1]])
         list(y=  eval(d[1]),
              m=eval(d[2]),
              d= eval(d[3]),
              h= eval(d[4]),
              m= eval(d[5]),
              s= eval(d[6]),
              tz=system.call("date '+%Z'"))
        }
     }
   if(is.S()) 
     {system.call <- function(cmd)  unix(cmd)            
      file.exists <- function(file) {1 == unix(paste("if [ -f ", file, 
           " ] ; then (echo 1) ; else (echo 0); fi"))}
      date.parsed <- function() 
        {d <- parse(text=system.call("date '+%Y %m %d %H %M %S'"),white=T)
         list(y=  eval(d[1]),
              m=eval(d[2]),
              d= eval(d[3]),
              h= eval(d[4]),
              m= eval(d[5]),
              s= eval(d[6]),
              tz=system.call("date '+%Z'"))
        }
     }
  }




##############################################################################

#    6/  Random number generation.

##############################################################################



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

        if(RNGkind()[2] == "default")  return(runif.default(n, min=min, max=max))
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

if (is.R())
  {
#  Prior to R 0.99 Wichmann-Hill was the default and the DSE version of
#  Box-Muller was used for rnorm.
  
    set.RNG <- function(kind=NULL, seed=NULL, normal.kind=NULL)
      {# with a null argument this also serves as get.RNG 
       #The next line means that set.RNG with null args does not truly
       #  return the state of the RNG in the case when it has not been 
       #  initialize. It first initializes it and then returns the state. The
       #  rational is that querying the state is usually for the purpose of
       #  reproducing it, so it must first be initialized to put it in a 
       #  reproducible state.
	if (!exists(".Random.seed")) z <- runif(1)
	old <- list(kind=RNGkind()[1],normal.kind=RNGkind()[2], seed=.Random.seed[-1])
        if (is.null(kind) & is.null(seed) & is.null(normal.kind)) return (old)
	if (is.list(kind)) 
          {seed        <- kind$seed
	   normal.kind <- kind$normal.kind
	   kind        <- kind$kind
	  }
	remove(".Random.seed", envir=.GlobalEnv) # otherwise RNGkind complains
	RNGkind(kind=kind, normal.kind=normal.kind)
	if ( 1==length(seed)) set.seed(seed) 
        else assign(".Random.seed", c(.Random.seed[1], seed), envir=.GlobalEnv)
	#RNGkind(kind=kind, normal.kind=normal.kind)
	old
      }

  }  # end of if is.R


#########################################################

#   test function

#########################################################


random.number.test <- function()
 {cat("Random number generator tests ...")
  if (is.R())  
     {test.seed<- 979   #previous to R 1.0.0: c( 979, 1479, 1542) 
      # values from 0.49 beta
      #test.valueU <-c(5.693354055333957e-01,1.051357751852140e-01,
      #    5.846933178718317e-02, 7.537960906527452e-02, 7.043734921992200e-01)
      #test.valueN <-c(-5.559389931781886e-01,
      #                   -1.902431069568611e+00,  1.524595894866778e+00,
      #                   -7.863494805034426e-01, 1.328128164898773e-01)
      # values from 0.99.0a
      # test.valueU <-c(0.25603057077527752, 0.07879165329010961,
      # 		 0.60394682330171257, 0.20843868707503158, 0.97636939375111098)
      # test.valueN <-c( -1.39965726956837644, -0.24025807684466990,
      #          2.34362137305187446, -0.66321208109989371, -0.71183559894654214)
      # values from 1.1.0 (devel) should also work in 1.0.1
      test.valueU <-c(0.59132864479704950, 0.76406894316060192,
                  0.18576870606880833, 0.81087542344137897, 0.05835627439859235)
      test.valueN <-c( 0.959409416509782953, 0.046546246156130192,
            -0.775306427558391964, -0.777761120325662803, -1.363043207314267313)
     }
  if (is.Splus()) 
     {test.seed<- c(37, 39, 39, 4, 7, 2, 27, 58, 38, 15, 32, 2)
      test.valueU <- c(0.4299328043125570, 0.3092006836086512,
            0.5808096211403608, 0.3061958812177181, 0.8137333435006440)
      test.valueN <- c( -0.7613318231781665, -0.5724360196433543,
            0.8536399448225964, -0.2269096022522968, -0.8126790170570223)
     }

  old.seed <- set.RNG(kind="default", seed=test.seed, normal.kind="default")
  on.exit(set.RNG(old.seed))

  ok <- TRUE
  if (1e-14 < max(abs(runif(5)-test.valueU)))
    {warning("The default runif number generator has been changed.")
     ok <- FALSE
    }

  set.RNG(kind="default", seed=test.seed, normal.kind="default")

  if (1e-14  < max(abs(rnorm(5)-test.valueN)))
    {warning("The default rnorm number generator has been changed.")
     ok <- FALSE
    }

  set.RNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
  if (set.RNG()$kind        != "Wichmann-Hill" ||
      set.RNG()$normal.kind != "Box-Muller"    ||
      all(set.RNG()$seed    != c(979,1479,1542) )) 
     {warning("RNG is not being set properly")
      ok <- FALSE
     }
  if (1e-14 < max(abs(runif(5) -
      c(0.56933540553339546, 0.10513577518521355, 0.05846933178718317,
        0.07537960906527452, 0.70437349219921996))))
    {warning("The Wichmann-Hill runif number generator has been changed.")
     ok <- FALSE
    }

# for the next R 1.0.0 
# set.RNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
# rnorm gives
#[1] -1.92425218107175877 -0.89568905204068128  2.12213361588187510
#[4]  0.81669202948845299 -0.13569189805256629
# as does
  set.RNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
  if (1e-14 < max(abs(rnorm(5) -
      c(-1.92425218107175877, -0.89568905204068128,  2.12213361588187510,
         0.81669202948845299, -0.13569189805256629)
      # pre R 1,0,0 c(0.4605069059114530, 0.7685565310963474, -0.3737680932387061,
      # 0.5926372779538560, 1.4995245125275518)
	)))
    {warning("The Box-Muller rnorm number generator has been changed.")
     ok <- FALSE
    }
  # this needs to be done twice with odd and even n to chech completely
  if (1e-14 < max(abs(rnorm(5) -
      c(-0.4602838255495997,  0.2761541652763959,  1.3265434523731297,
        0.6856247181400722, -1.8336523890846541) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }
  if (1e-14 < max(abs(rnorm(6) -
      c(1.9850437759531543,  0.6107700961454843, -0.9419893721776470,
        1.1031328847642050,  0.4184702210057414,  0.9167797157851526) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }

  if (1e-14 < max(abs(rnorm(6) -
      c(-0.724539745179790251, -0.439138566092752758,  1.466237618877826110,
         0.289289597398559639,  0.003007778996985022,  1.008712871048744297) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }

  if (ok) cat("ok\n")
  else    cat("failed!\n")
  invisible(ok)
 }



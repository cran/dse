##########

# Functions for S should be loaded from the S/ directory before 
#    functions in this file.

###########


#  Functions are for identifying S or R and flavours.
  
#  Functions depending on differences between S and R
	 
#  Random number generation.
#    Done with the aid of example from B. D. Ripley.

#if (!exists("is.R"))    is.R <- function() {FALSE}

# This is only needed for "For" loops which are never called in R versions
#if(is.R()) synchronize <- function(){NULL} # perhaps this should do something?

# For S different versions override later 
# if(is.R()) 
   {is.S <- is.Svanilla <- is.Splus <- is.Splus.pre3.3 <- function(){FALSE}

    file.date.info <- function(file)
     	  {# format of the result could be better. 
	   x <- as.POSIXlt(file.info(file)$mtime)
	   c(1+x$mon,x$mday,x$hour,x$sec) # as previously used. should be improved
     	  }
   
    date.parsed <- function() 
          {d <- as.POSIXlt(Sys.time())
           list(y = 1900 + d$year,   # as previously used. should be improved
              m = 1+ d$mon,
              d = d$mday,
              H = d$hour,
              M = d$min,
              S = d$sec,
              tz = attr(d, "tzone"))
          }
        #  syskern.rm() was previously unlink() which is now supported in R.
        #  Unfortunately the argument recursive is not supported in S and the
        #    R 1.2 default value of FALSE is a change from previous Unix versions of R
        #    and from S and causes problems the way it is used in DSE.
    syskern.rm <- function(file) unlink(file, recursive = TRUE)
	
    .SPAWN <- FALSE

 }

if (!exists("Sys.mail")) { 

Sys.mail <- function(address = Sys.info()$user,
                     ccaddress  = NULL,
                     bccaddress = NULL,
		     subject    = NULL,
		     body= "no body",
                     method = getOption("mailer")) {
   if(!(is.character(address) && nchar(address)>0))
      stop("A character string address must be specified.")

   # The arguments to mail, mailx, and Mail are all the same, but a different
   # mailer will require that this first part be re-organized under the
   # specific method.
   file <- tempfile()
   on.exit(unlink(file))
   cat(body, file=file, append=FALSE, sep="\n")
   cmdargs <- paste(address, "<", file, "2>/dev/null")
	
   if(is.character(ccaddress) && nchar(ccaddress)>0) 
            cmdargs <- paste(" -c '", ccaddress, "' ",  cmdargs)

   if(is.character(bccaddress) && nchar(bccaddress)>0) 
            cmdargs <- paste(" -b '", bccaddress, "' ",  cmdargs)

   if(is.character(subject) && nchar(subject)>0) 
            cmdargs <- paste(" -s '", subject, "' ",  cmdargs)

   status <- 1
   if(method == "mailx") status <- system(paste("mailx", cmdargs)) else
   if(method == "mail") status <- system(paste("mail", cmdargs))   else 
   if(method == "Mail") status <- system(paste("Mail", cmdargs))   else {
	warning(paste("mail method ", method, " not supported.\n"))
	return(FALSE)
	}
   if(status > 0) {
	     warning(paste("sending email with ", method, " failed.\n"))
	     return(FALSE)
	     }
   TRUE
   }

}
   
###########################################################

#    Random number generation.

###########################################################



 
#if (is.R())  In S different versions overrid later
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

#   test function (This test is run by other tests so it should
#   be defined here and not moved to tests directory.

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

  if (ok) {cat("ok\n"); invisible(TRUE) } else { cat("failed!\n"); stop("failed")}
 }
  

\name{00General_Description}
\alias{00General_Description}
\title{Syskern}
\center{OS and S/R independent programming kernel}
\usage{
	library("syskern")}
\description{Kernel of functions for OS and S/R independence, including RNG. ...}
\details{
This library has code which contains operating system and S/R specific 
functions. They are intended to be used as a kernel to help 
protect other code from these problems. In many cases the functionality is 
supplied by other functions in S or R, and in some cases the 
function names have been used in R. The purpose of this package is to provide a 
common interface so that other code can be written to work in S, R, and on 
different operating systems.

MSwindows versions are not done. If
anyone can contribute these please contact <pgilbert@bank-banque-canada.ca>.

S/R differences are handled by functions global.assign, system.info, 
exists.graphics.device, synchronize, and list.add for [["new.element"]]<-.
      
The library also implements an approach to random number generation 
which allows the same random experiments to be replicated in S and R.
As of R 1.0.1 this adds nothing to R except the function set.RNG which sets
both RNGkind and the seed. The functions allow the R results using 
Wichmann-Hill and Box-Muller to be replicated in S.
These functions are intended primarily as a way to confirm that other code works in 
the same way in both S and R, not as an improved RNG. Default and other
RNGs can  still be used and are probably preferred for both speed and 
theoretical reasons. The present version of this works only for uniform 
and normal distributions.

OS differences are handled by functions system.call, sleep, 
present.working.directory, whoami, file.copy, 
file.date.info, date.parsed, mail, local.host.netname, 

Also a number of is.XXX functions are defined to identify systems.

 The variable  .SPAWN is also set to be used to identify if Splus "For" loops
    should be used. (It is sometimes better not to use these even in Splus.)
}
/section{General Logic and organization of these functions} {

 1/ The first group of functions are for identifying S or R and flavours.
 
 2/ The second group of functions are for identifying the operating system. 
 
 3/ The third group specify a few functions which depend only on the 
         differences between S and R.
	 
 4/ The fourth group specify functions which depend only on the 
         differences among operating system.
	 
 5/ The fifth group specify a few functions which depend on both R/S and the 
         differences among operating system.
	 
 6/ The sixth group are functions enable the same random number generation
         (uniform and normal) in Splus and R. This was done with the aid of 
	 some example from B. D. Ripley.

  >>> I would very much like any input WRT  MS Windows / Win95 / NT / Mac <<<
   
 The function system.call is defined in order to provide a generic way to
  make a call to the operating system. When the calls are specific
  to Unix then the function unix() might be used (though that is now 
  deprecated in R and produces a warning messsage). However, in general these
  kernel functions are not to give a generic way to make any call to   
  the operating system, but rather a generic way to make some kinds of specific  
  requests to the operating system (like date, mail, sleep, whoami).
}

\examples{
  is.R()
  is.unix()
  set.RNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
  rnorm(10)
}
\keyword{programming}
\keyword{interface}
\keyword{utilities}
\keyword{distribution}

\author{Paul Gilbert}


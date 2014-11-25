\name{00General.Description}
\alias{00General.Description}
\title{Syskern}
\center{OS and S/R independent programming kernel}
\usage{
	library("syskern")}
\description{Kernel of functions for OS and S/R independence, including RNG. ...}
\details{
The purpose of this package is to provide a 
common interface so that other code can be written to work in S, R, and on 
different operating systems.
The package has functions which provide wrappers for operating system requests
and for functions which work differently in various versions of S and/or R. 
The wrappers are intended to be used as a kernel of functions to help 
protect other code from these problems. In many cases the functionality is 
supplied by other functions in S or R, and in some cases the 
function names have been used in R.
Indeed, as of R 1.2 many of the functions
are longer be necessary in R and the libraries main purpose is to provide
similar functionality in S. 

MSwindows versions are not working for S, but as noted above, R provides many of the
functions for all operating systems. If anyone can contribute to the S version
then these please contact \email{pgilbert@bank-banque-canada.ca}.

S/R differences are handled by functions \code{global.assign},
\code{syskern.rm} (for \code{unlink}),
\code{system.info}, \code{exists.graphics.device}, and \code{list.add} 
for \code{[["new.element"]]<-}.

The functions \code{dev.ask} and \code{synchronize}(does nothing) are defined 
for R, and for S the functions \code{Sys.getenv}, \code{getwd}, and
\code{require} (does nothing) are defined.
      
The library also implements an approach to random number generation 
which allows the same random experiments to be replicated in S and R.
As of R 1.0.1 this adds nothing to R except for the function set.RNG which sets
both RNGkind and the seed. The functions allow the R results using 
Wichmann-Hill and Box-Muller to be replicated in S.
These functions are intended primarily as a way to confirm that other code works in 
the same way in both S and R, not as an improved RNG. Default and other
RNGs can  still be used and are probably preferred for both speed and 
theoretical reasons. The present version of this works only for uniform 
and normal distributions. These were done with the aid of example 
from B. D. Ripley.


OS differences are handled by functions \code{system.call}, \code{Sys.sleep}, 
\code{getwd}, \code{whoami}, \code{file.copy}, 
\code{file.date.info}, \code{date.parsed}, \code{mail} and 
\code{local.host.netname}. 

Also a number of \code{is.XXX} functions are defined to identify systems.

The variable \code{.SPAWN} is also set to be used to identify if Splus "For" loops
    should be used. (It is sometimes better not to use these even in Splus.)
}
  
 The function system.call is defined in order to provide a generic way to
  make a call to the operating system. When the calls are specific
  to Unix then the function unix() might be used (though that is now 
  deprecated in R and produces a warning messsage). However, in general these
  kernel functions are not to give a generic way to make any call to   
  the operating system, but rather a generic way to make some kinds of specific  
  requests to the operating system (like date, mail, Sys.sleep, whoami).


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


\name{DSEflags}
\alias{.DSEflags}
\title{Flags to Indicate Use of Compiled Code}
\description{
    Determines if compiled code should be used or not.
    }
\usage{
    .DSEflags(new)
    .DSEflags()
    }
\arguments{
    \item{new}{A list which must have elements COMPILED and DUP.}
     }
\sideffects{
    Setting flags with this function is primarily for debugging. It should not
    normally be needed by users.
    Several DSE functions which call compiled fortran or C code will use
    the equivalent S/R version if the .DSEflags()$COMPILED returns 
    FALSE.  .DSEflags()$DUP is passed to
    fortran calls and controls whether R duplicates arguments passed
    to the fortran code. The safe setting is TRUE. Setting FALSE saves some
    memory but does not seem to give much speed gain. 
    }
\examples{
    .DSEflags(list(COMPILED=TRUE, DUP=TRUE))
    .DSEflags()$COMPILED
}
\concept{DSE}
\keyword{programming}



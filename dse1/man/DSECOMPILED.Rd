\name{DSECOMPILED}
\alias{.DSECOMPILED}
\alias{.DSEDUP}
\title{Use Compiled or S/R Code}
\description{
    Determines if compiled code should be used or not.
    }
\usage{
    .DSECOMPILED <- TRUE
    .DSECOMPILED <- FALSE
    .DSEDUP <- TRUE
    .DSEDUP <- FALSE
}

\sideffects{
    Several DSE functions which call compiled fortran or C code will use
    the equivalent S/R version if the .DSECOMPILED global variable is set 
    to FALSE. Its main purpose is for debugging. .DSEDUP is passed to
    fortran calls and in R controls whether R duplicates arguments passed
    to the fortran code. The safe setting is TRUE. Setting FALSE saves some
    memory but does not seem to give much speed gain. 
}
\examples{
    .DSECOMPILED <- TRUE
    .DSEDUP <- TRUE
}
\concept{DSE}
\keyword{programming}



\name{freeze}
\alias{freeze}
\alias{freeze.default}
\title{Get fixed data snapshot}
\description{Generic method. See specific methods for details.}

\usage{
    freeze(data, ...)
    \method{freeze}{default}(data, ...)
}

\arguments{
    \item{data}{A data source description.}
    \item{...}{Additional arguments for specific methods.}
}

\details{
This function extracts data from 
a database (for example using the TS PADI programs are available 
at \url{www.bank-banque-canada.ca/pgilbert}). This method is generic.
The function logically belongs with the tfpadi functions in the dsepadi package
but the generic and default are included here since they are required in many
functions whereas the database interface has much more limited use.
Typically the argument data is a tfPADIdata or TSPADIdata object identifying 
the source of the data. See help for tfPADIdata and TSPADIdata. 
The default method usually just
returns its argument, so freeze has no effect. This way freeze
can be used to write functions which will take a snapshot from the database
when they execute or will work with an already fixed copy of data if that is
what is supplied. The default does allow for a character argument, in which case
it is used to construct a tfPADIdata object using server="ets", then freeze that
object. This allows for a simple syntax to grab a series from the database. The
server="ets" is for convenience at the Bank of Canada and another default server
might be more convenient elsewhere.
}

\seealso{ 
    \code{\link[dsepadi]{tfPADIdata}}
    \code{\link[dsepadi]{tfputpadi}}
    \code{\link[dse2]{freeze.tfPADIdata}}
    \code{\link[dsepadi]{freeze.TSPADIdata}}
    \code{\link[dse2]{TSPADIdata}}
}
\keyword{programming}
\keyword{ts}
\keyword{data}



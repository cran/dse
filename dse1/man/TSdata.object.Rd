\name{TSdata.object}
\alias{TSdata.object}
\title{time series data object}
\description{
Class "TSdata" of time series data objects for use with TSmodels.
}

\section{Generation}{
This class of objects is returned by specific methods of the 
function TSdata or can be built according to the description below.
}

\section{Methods}{
The \code{TSdata} class of objects has methods for the generic
functions \code{print, plot, periods, start, end, ..., test.equal, series.names,
input.series.names, output.series.names}
Also, the function \code{is.TSdata} is supported.
}

\section{Inheritance}{
Other data classes inherit from the class \code{TSdata}.}

\section{Structure}{
Objects are a list with class the most general class "TSdata". The
native form for this library has elements \code{input} and
\code{output}. Any other elements are ignored. \code{input} and
\code{output} are matrices (or tframe or time series matrices) of the
input and output data, with each series in a column.  \code{TSPADIdata}
objects inherit from this class but have a somewhat different
structure.  \code{TSPADIdata} makes it possible to retrieve data from
an external database when it is needed. These subclass objects do not
contain the actual data, but only the names of the series and the data
base where they are located. The function set.TSPADIdata can be used to set
up an object of class c("TSPADIdata", TSdata").
}

\seealso{
\code{\link{TSdata}}
\code{\link{TSmodel}}
\code{\link{TSestModel.object}}
\code{\link[dsepadi]{TSPADIdata}}
}
\keyword{DSE}
\keyword{ts}


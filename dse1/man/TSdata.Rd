\name{TSdata}
\alias{TSdata}
\alias{TSdata.default}
\alias{TSdata.TSdata}
\alias{TSdata.TSestModel}
\alias{is.TSdata}
\alias{as.TSdata}

\title{Construct TSdata time series object}
\description{
    Constructor for generating or extracting a "TSdata" object 
    containing data for use by TSmodels.
}
\usage{
    TSdata(data, ...)
    \method{TSdata}{default}(data, input=NULL, output=NULL, ...)
    \method{TSdata}{TSdata}(data, ...)
    \method{TSdata}{TSestModel}(data, ...)
    is.TSdata(obj)
    as.TSdata(d)
}
\details{
    Generic method to construct or extract a TSdata object. The default 
    method constructs a TSdata object. Specific methods extract the TSdata
    from other objects (which must contain TSdata).
    The function \code{is.TSdata(data)} returns TRUE is data inherits 
    from "TSdata" and FALSE otherwise.

    The function \code{as.TSdata} uses the elements input and output 
    directly and
    strips away other class information and parts of the object (and does 
    not make use of input.data(data) or output.data(data) which may 
    do something special for certain classes.
}
\seealso{
\code{\link{TSdata.object}},
\code{\link[tframe]{freeze}},
\code{\link[dse2]{freeze.TSPADIdata}},
\code{\link{TSmodel}},
\code{\link{TSestModel.object}}
\code{\link[dse2]{set.TSPADIdata}}
}
%\keyword{DSE}
\keyword{ts}


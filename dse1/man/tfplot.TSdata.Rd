\name{tfplot.TSdata}
\alias{tfplot.TSdata}
\alias{tfplot.TSmodel}
\alias{tfplot.TSestModel}

\title{Tfplot Specific Methods}
\description{See the generic function description.}
\usage{
    \method{tfplot}{TSdata}(..., start.=NULL,end.=NULL, Title="", 
        reset.screen=T,
        select.inputs =seq(length= input.dimension(data)),
        select.outputs=seq(length=output.dimension(data)),
	mar=par()$mar, 
        graphs.per.page=5, ylab=NULL)
    \method{tfplot}{TSestModel}(..., start.=NULL,end.=NULL,Title=NULL, 
        reset.screen=T, select.inputs=NULL, select.outputs=NULL,
        graphs.per.page=5, mar=par()$mar)
    }

\seealso{ \code{\link{tfplot}} }

%\keyword{DSE}
\keyword{ts}


\name{test.Equaltframes}
\alias{test.Equaltframes}
\alias{test.Equaltframes.default}
\alias{test.Equaltframes.stamped}

\title{Compare Two Time Frames}
\description{
Generic function to compare two time frames.  The methods return a 
logical value, TRUE if the time frames are the same type and value and FALSE otherwise.
}
\usage{
    test.Equaltframes(tf1,tf2)
    \method{test.Equaltframes}{default}(tf1,tf2)
    \method{test.Equaltframes}{stamped}(tf1,tf2)
    }
\arguments{
    \item{tf1, tf2}{Time frames of the same class.}
}
\value{TRUE or FALSE}
\details{
Time frames are compared. Time frames need to be of the same class (although it would be nice if they did
not need to be).
}
\examples{
test.Equaltframes(tframe(matrix(1:10,10,2)), tframe(array(1:10, c(10,2))))
}
\seealso{
    \code{\link{test.equal}}
}
\keyword{programming}
\keyword{utilities}
\keyword{ts}
\keyword{chron}



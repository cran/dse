\name{syskern.rm}
\alias{syskern.rm}
\alias{file.date.info}
\alias{date.parsed}

\title{Utilities to handle small S/R differences}
\description{Utilities to handle small S/R differences.}

\usage{
    date.parsed()
    syskern.rm(file)
    file.date.info(file)
}

\arguments{
    \item{file}{a character string indicating a file name.}
}
\value{misc.}
\details{
    Utilities to handle small S/R differences. 
    \code{syskern.rm} does unlink with default recursive=TRUE.
}
\seealso{
    \code{\link{set.RNG}}
    \code{\link{is.S}}
    \code{\link{file.copy}}
    \code{\link{Sys.info}}
}
\examples{date.parsed() }
\keyword{programming}
\keyword{utilities}



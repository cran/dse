\name{syskern.rm}
\alias{syskern.rm}
\alias{file.date.info}
\alias{date.parsed}
\alias{synchronize}

\title{Utilities to handle small S/R differences}
\description{Utilities to handle small S/R differences.}

\usage{
    date.parsed()
    synchronize()
    syskern.rm(file)
    file.date.info(file)
    Sys.getenv()
}

\arguments{misc.}
\value{misc.}
\details{
    Utilities to handle small S/R differences. 
    \code{syskern.rm} does unlink with default recursive=T.
    \code{synchronize} does nothing in R.
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



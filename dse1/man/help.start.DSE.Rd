\name{help.start.DSE}
\alias{help.start.DSE}
\title{Start DSE Help in S}
\description{Start a browser and load DSE help.}
\usage{
    help.start.DSE(browser = "netscape")
}
\arguments{
    \item{browser}{A system command for starting a browser.}
}
\value{None.}
\details{
    The DSE HTML help start page is loaded by the specified browser. 
    This function is for use in S. In R the help is integrated with R's
    help facilities.
}
\sideffects{A browser is started and DSE HTML help is loaded.}
\examples{if(is.S()) help.start.DSE()}
%\keyword{DSE}
\keyword{ts}



\name{summary}
\alias{summary}
\title{Display Summary Information.}
\description{
Summary is a generic  function which can be used on DSE objects.
}
\usage{
	summary(object)
}
\arguments{
	\item{obj}{An object of class TSdata, TSmodel or TSestModel.}
}
\value{No value is returned}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	model <- est.VARX.ls(eg1.DSE.data.diff)
	summary(model)
}
\keyword{DSE}
\keyword{ts}


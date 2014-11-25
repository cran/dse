\name{scale.TSdata}
\alias{scale.TSdata}
\alias{scale.TSmodel}
\alias{scale.TSestModel}
\alias{scale.default}
\alias{scale.ARMA}
\alias{scale.innov}
\alias{scale.non.innov}
\alias{scale.check}
\alias{scale.check.TSestModel}
\alias{scale.check.TSmodel}

\title{Scale Methods for TS objects}
\description{Scale data or a model by a given factor.}
\usage{
	scale(obj, scale)}
\arguments{
	\item{obj}{TSdata, TSmodel or an object containing these.}
	\item{scale}{
	A list with two matrices or vectors, named input and output,
	giving the multiplication factor for inputs and outputs.
	Vectors are treated as diagonal matrices.
	\code{scale$input} can be NULL if no transformation is to be applied (or 
	the data or model has no input.)}
}

\value{
	The resulting data or model is different from
	the original in proportion to scale.  ie. if S and T are output and 
	input scaling matrices then 
         y'(t) = S y(t) where y' is the new output
         u'(t) = S u(t) where u' is the new input

	For models the result has inputs and outputs (and innovations) 
	which are scaled as if data scaling had been applied to them as above. 
	Thus if the input and output scales are diagonal matrices or 
	scalars the plot of the predictions and residuals for
	l(scale(model,somescale), scale(data, somescale))
	while have the same appearance as l(model, data) but will be 
	scaled differently.
}
\seealso{
	\code{\link[base]{scale}}
}
\examples{
	if(is.R()) data("eg1.DSE.data.diff", package="dse1")
	# This is a simple example. Usually scale would have something 
	# to do with the magnitude of the data.
	z <- scale(eg1.DSE.data.diff, 
	    list(input=rep(2, input.dimension(eg1.DSE.data.diff)), 
		output=rep(2,output.dimension(eg1.DSE.data.diff))))
	model <- est.VARX.ls(eg1.DSE.data.diff)
	model <- scale(model, 
	    list(input=rep(2, input.dimension(eg1.DSE.data.diff)), 
		output=rep(2,output.dimension(eg1.DSE.data.diff))))
}
\keyword{DSE}
\keyword{ts}


\name{min.forecast.cov}
\alias{min.forecast.cov}
\title{Minimum forecast cov models}
\description{Extract the minimum forecast.cov at each horizon}
\usage{
	min.forecast.cov(obj, select.series=1, verbose=T)
}
\arguments{
\item{ obj}{An object as returned by mine.strip.}
\item{ select.series}
{An indicator of the series which are to be used as the bases for selection.}
\item{ verbose}{If true additional information is printed.}
}
\value{
The returned object is a vector indicating the element of forecast.cov which
was the min at each horizon.
}
\details{
Select the min covariance (for select.series only!) at each horizon and
print.  The returned object is a vector indicating the element of
forecast.cov which was the min at each horizon. It is suitable as an
argument to plot eg:  tfplot(obj, select.cov=min.forecast.cov(obj)) The
results of this plot are similar to the default results of
tfplot(select.forecast.cov(obj)). Covariance information and information
about the horizon where the model is optimal are given.
}
\seealso{
\code{\link{select.forecast.cov}},
\code{\link{exclude.forecast.cov}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <- mine.strip(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(est.VARX.ls=list(max.lag=3)))
z <-  min.forecast.cov(z)
}
\keyword{DSE}
\keyword{ts}


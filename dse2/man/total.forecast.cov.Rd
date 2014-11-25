\name{total.forecast.cov}
\alias{total.forecast.cov}
\title{Sum covariance of forecasts across all series}
\description{Sum covariance of forecasts across all series.}
\usage{
    total.forecast.cov(obj, select=NULL)
}
\arguments{
    \item{obj}{An object  as returned by forecast.cov.}
    \item{select}{
    Series to be select for summation. With the default all series are selected.}
}
\value{
An object similar to that returned by forecast.cov, with the covariance
summed over all selected series.
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
model1 <- est.VARX.ar(eg1.DSE.data.diff)
model2 <- est.VARX.ls(eg1.DSE.data.diff)
z <-  total.forecast.cov(forecast.cov(model1, model2,
                         data=trim.na(eg1.DSE.data.diff)))
}
%\keyword{DSE}
\keyword{ts}


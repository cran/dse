\name{est.black.box3}
\alias{est.black.box3}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    est.black.box3(data, estimation='est.VARX.ls', 
          lag.weight=1.0, 
          reduction='reduction.Mittnik', 
          criterion='aic', 
          trend=F, 
          subtract.means=F,  re.add.means=T, 
          standardize=F, verbose=T, max.lag=12, sample.start=10)
}
\arguments{
    \item{data}{A TSdata object.}
    \item{estimation}{A character string indicating the estimation method to use.}
    \item{lag.weight}{Weighting to apply to lagged observations. }
    \item{reduction}{Character string indicating reduction procedure to use. }
    \item{criterion}{
    Character string indicating model selection criteria. taic might be a
    better default selection criteria but it is not available for ARMA
    models.}
    \item{trend}{If T include a trend in the model. }
    \item{subtract.means}{
    If T the mean is subtracted from the data before estimation.}
    \item{re.add.means}{
    If subtract.means is T then if re.add.means is T the estimated model is
converted back to a model for data without the mean subtracted.}
    \item{standardize}{
    If T the data is transformed so that all variables have the same variance.}
    \item{verbose}{
    If T then additional information from the estimation and reduction 
procedures is printed.}
    \item{max.lag}{The number of lags to include in the VAR estimation.}
    \item{sample.start}{The starting point to use for calculating information criteria.}
}
\value{A TSestModel.}
\details{
VAR models are estimated for each lag up to the specified max.lag. From 
these the best is selected according to the specified criteria. The 
reduction procedure is then applied to this best model and the best 
reduced model selected. The default estimation procedure is least 
squares estimation of a VAR model.
}
\seealso{
\code{\link{est.black.box1}},
\code{\link{est.black.box2}}
\code{\link{est.black.box4}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  est.black.box3(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


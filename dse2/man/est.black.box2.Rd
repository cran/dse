\name{est.black.box2}
\alias{est.black.box2}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    est.black.box2(data, estimation='est.VARX.ls', 
          lag.weight=.9, 
          reduction='reduction.Mittnik', 
          criterion='taic', 
          trend=F, 
          subtract.means=F,  re.add.means=T, 
          standardize=F, verbose=T, max.lag=12)
}
\arguments{
    \item{data}{A TSdata object.}
    \item{estimation}{A character string indicating the estimation method to use.}
    \item{lag.weight}{Weighting to apply to lagged observations. }
    \item{reduction}{Character string indicating reduction procedure to use. }
    \item{criterion}{Character string indicating model selection criteria. }
    \item{trend}{If T include a trend in the model. }
    \item{subtract.means}{
    If T the mean is subtracted from the data before estimation.}
    \item{re.add.means}{
    If subtract.means is T then if re.add.means is T 
the estimated model is
converted back to a model for data without the mean subtracted.}
    \item{standardize}{
    If T the data is transformed so that all variables have the same variance.}
    \item{verbose}{If T then additional information from the estimation and 
reduction procedures is printed.}
    \item{max.lag}{The number of lags to include in the VAR estimation.}
}
\value{A TSestModel.}
\details{
A model is estimated and then a reduction procedure applied. The
default estimation procedure is least squares estimation of
a VAR model with lagged values weighted. This procedure is discussed in
Gilbert (1995).
}
\references{
    Gilbert, P.D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.
}
\seealso{
\code{\link{est.black.box1}},
\code{\link{est.black.box3}}
\code{\link{est.black.box4}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  est.black.box2(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


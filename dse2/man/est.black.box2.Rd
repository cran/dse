\name{est.black.box2}
\alias{est.black.box2}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel.}
\usage{
    est.black.box2(data, estimation='est.VARX.ls', 
          lag.weight=.9, 
          reduction='reduction.Mittnik', 
          criterion='taic', 
          trend=FALSE, 
          subtract.means=FALSE,  re.add.means=TRUE, 
          standardize=FALSE, verbose=TRUE, max.lag=12)
}
\arguments{
    \item{data}{a TSdata object.}
    \item{estimation}{a character string indicating the estimation method to use.}
    \item{lag.weight}{weighting to apply to lagged observations. }
    \item{reduction}{character string indicating reduction procedure to use. }
    \item{criterion}{criterion to be used for model 
       selection. see \code{information.tests.calculations}.}
    \item{trend}{if TRUE include a trend in the model. }
    \item{subtract.means}{
    if TRUE the mean is subtracted from the data before estimation.}
    \item{re.add.means}{
    if subtract.means is TRUE then if re.add.means is TRUE 
the estimated model is
converted back to a model for data without the mean subtracted.}
    \item{standardize}{
    if TRUE the data is transformed so that all variables have the same variance.}
    \item{verbose}{if TRUE then additional information from the estimation and 
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
\code{\link[dse1]{information.tests.calculations}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  est.black.box2(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}


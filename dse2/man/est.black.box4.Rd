\name{est.black.box4}
\alias{est.black.box4}
\alias{bft}
\title{Estimate a TSmodel}
\description{Estimate a TSmodel with Brute Force Technique.}
\usage{
    est.black.box4(data, estimation="est.VARX.ls", 
                lag.weight=1.0,  variable.weights=1, 
                reduction="reduction.Mittnik", 
                criterion="taic", 
                trend=F, subtract.means=F,  re.add.means=T, 
                standardize=F, verbose=T, max.lag=12, sample.start=10, warn=T)
    bft(data, ... )
}
\arguments{
    \item{data}{A TSdata object.}
    \item{estimation}{A character string indicating the estimation method to use.}
    \item{lag.weight}{Weighting to apply to lagged observations. }
    \item{variable.weights}{
    Weighting to apply to series if estimation method is est.wt.variables. }
    \item{reduction}{Character string indicating reduction procedure to use. }
    \item{criterion}{Character string indicating model selection criteria. }
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
    \item{max.lag}{VAR estimation is done for each lag up to max.lag.}
    \item{sample.start}{
    The starting point to use for calculating information criteria in the 
    final selection.}
}
\value{A TSestModel.}
\details{
For each lag up to max.lag a VAR model is estimated and then a
reduction procedure applied to select the best reduced model. Finally
the best of the best reduced models is selected. The default estimation
procedure is least squares estimation of the VAR models. This procedure
is described as the brute force technique (bft) in Gilbert (1995).
}
\references{
    Gilbert, P.D. (1995) "Combining VAR Estimation and State Space 
    Model Reduction for Simple Good Predictions" J. of Forecasting: 
    Special Issue on VAR Modelling. 14:229-250.
}
\seealso{
\code{\link{est.black.box1}},
\code{\link{est.black.box2}}
\code{\link{est.black.box3}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  bft(eg1.DSE.data.diff)
}
%\keyword{DSE}
\keyword{ts}



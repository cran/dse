\name{00General.Description}
\alias{00General.Description}
\title{DSE}
\center{Dynamic Systems Estimation - Multivariate Time Series Package}
\usage{
   library("dse1") # 1 is the main part
   library("dse2") #   extensions
   library("dse3") #   extensions
  }

\description{
This is a library of functions for time series modeling.  The library 
works with the S (Splus) and R languages. A "Brief User's Guide" is available 
at \url{www.bank-banque-canada.ca/pgilbert}  
This help and other information is also available at that web location.

The library implements an object oriented approach to time series
modeling. This means that different
model and data representations can be implemented with fairly simple 
extensions to the library. 

The library includes multi-variate state 
space and ARMA (including VAR) models. The library also implements 
Troll models as an example of another class of model. (Troll models are run by completely 
separate software, Troll from Intex Inc, as if they were an integral part 
of the library. Models and data are passed back and forth to Troll and the
results can be analyzed with already existing functions in the
library. Note: the Troll interface is broken in the current version.)

The library includes methods for simulating, estimating, and converting 
among different
model representations. These form the basic part of the library. Methods for
studying estimation techniques and for examining the forecasting properties
of models are also documented in the User's Guide.

There are also functions for forecasting and for
evaluating the performance of forecasting models as well as 
functions for evaluating model estimation techniques.

DSE requires libraries tframe and syskern. The first provides a kernel of 
functions for 
programming time series methods relatively independently of the representation 
of time. The second provides some utilities to handle S/R differences 
and OS differences. For programmers, these libraries may of interest 
by themselves. Tframe is intended to make it easier to write code which can use any new/better 
time representations when they appear. It also provides plotting, time 
windowing, and some other utility functions which are specifically intended 
for time series.

Relative to commercial packages the library is especially useful for
time series research (such as studying estimation methods). 
For usual time series applications there may be commercial
packages which are preferable.

The functions described in the
Brief User's Guide should work fairly reliably, however, many of the functions
described in the help facility are still under development and may not
work. In addition, there may be functions described in the help
facility for which the code is not yet distributed. This is a compromise which 
allows me to make the software
available with minimum effort. This software is not a commercial
product. It is the by-product of an ongoing research effort.
Constructive suggestions and comments are welcomed. I can be reached at
\email{pgilbert@bank-banque-canada.ca} or by phone at (613) 782-7346.

There is also a mailing list available which I hope will let users help 
one another, as I cannot always be available. Please subscibe if you use
the library extensively, especially if you can help other users.
To subscribe send mail to 
\email{boc\_list@bank-banque-canada.ca}
with the body line \code{subscribe boc\_dse Your Name}
}

\section{MainDataObjects}{
The main objects are classes "TSdata", "TSmodel", and  "TSestModel". 
These are, respectively,  time series data objects, models, and objects with 
a model, data and some estimation information. For each of these there are 
several sub classes.
}

\seealso{
\code{\link{TSdata}},
\code{\link{TSmodel}},
\code{\link{TSestModel.object}}
}
\keyword{DSE}
\keyword{ts}


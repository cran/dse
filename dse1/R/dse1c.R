#   2000/04/20 14:50:53
# For installation instructions see the file read.me or the brief user's
#    guide (postscipt file guide.ps).

##############################################################################



tfwindow.TSdata <- function(x, start=NULL, end=NULL, warn=T)
{# window a TSdata object
  if (0 != input.dimension(x))  input.data(x) <- 
      tfwindow(input.data(x), start=start, end=end, warn=warn)
  if (0 != output.dimension(x)) output.data(x) <- 
      tfwindow(output.data(x), start=start, end=end, warn=warn)
  x
}

window.TSdata <- tfwindow.TSdata 


combine <- function(e1,e2)UseMethod("combine")
combine.default <- function(e1,e2){list(e1,e2)}


combine.TSdata <- function(e1,e2)
{# make a new TSdata object with the two objects
 output.data(e1) <- tbind(output.data(e1),output.data(e2)) 
 if ((0 != (input.dimension(e1))) & (0 != (input.dimension(e2))))
       input.data(e1) <- tbind(input.data(e1),input.data(e2)) 
 else  {if (0 != (input.dimension(e2)))  input.data(e1) <-input.data(e2)  }
 e1
}


trim.na.TSdata <- function(data, start.=T, end.=T)
{# trim NAs from the ends of TSdata.
 # (Observations for all series are dropped if any one contains an NA.)
 # if start.=F then beginning NAs are not trimmed.
 # If end.=F   then ending NAs are not trimmed.
 # The same truncation is applied to both input and output
 p <- output.dimension(data)
 m <- input.dimension(data)
 if (m==0)
   mat <- trim.na(output.data(data))
 else
   mat <- trim.na(tbind(input.data(data),output.data(data)),start.=start.,end.=end.)
 tf <- tframe(mat)
 if (m!=0)
   {sn <- input.series.names(data)
    input.data(data)  <- tframed(mat[,1:m, drop=F], tf) 
    input.series.names(data) <- sn
   }
 sn <- output.series.names(data)
 output.data(data) <- tframed(mat[,(m+1):(m+p), drop=F], tf)
 output.series.names(data) <- sn
 data
}



diff.log <- function(x,  lag = 1, base = 2.71828182845905)
{#Calculate the difference from lag periods prior for log of data.
diff(log(x, base =base), lag=lag)
}


ytoypc <- function(ser) {
  # Convert level data to year over year percent change.
  # note: percent.change can alter the name, so grab it first.
  nm <- paste("y to y %ch", series.names(ser))
  ser <- percent.change(ser, lag=frequency(ser))
  series.names(ser) <- nm
  ser
 }
 


percent.change <- function(obj, ...) UseMethod("percent.change")

percent.change.list <- function(..., base=NULL, lag=1, cumulate=F, e=F)
  {#Calculate the percent change relative to the data lag periods prior.
   #... should be a list of objects to which percent.change can be applied.
   pchange <- list()
   for (mat in list(...))
          pchange <- append(pchange,list(percent.change(mat)))
   pchange
  }

percent.change.default <- function(mat, base=NULL, lag=1, cumulate=F, e=F)
{#Calculate the percent change relative to the data lag periods prior.
 # mat should be  a  matrix or vector.
 # If cumulate is T then the data is cumulated first. cumulate can be
 # a logical vector with elements corresponding to columns of m.
 # If e is T the exponent of the series is used (after cumulating 
 #   if cumulate is T). e can be
 # a logical vector with elements corresponding to columns of m.
 # If base is provided it is treated as the first period value 
 #  (prior to differencing). It is prefixed to the m prior to 
 #  cumulating. It should be a vector of length dim(m)[2]. 
 #  (If e is T then base should be log of the original data).
   cls <- tfclass(mat)
   if (is.tframed(mat)) tf <- list(end=end(mat), frequency=frequency(mat))
   else tf <- NULL
   if (is.null(dim(mat)))
     {vec <- T
      mat <- matrix(mat, length(mat),1)
     }
   else vec <- F
   mm <- rbind(base,mat)
   if (any(cumulate))
          mm[,cumulate] <-apply(mm[,cumulate,drop=F],2,cumsum)
   if (any(e)) mm[,e] <- exp(mm[,e,drop=F])
   N <- dim(mm)[1]
   pchange <-100*(mm[(lag+1):N,,drop=F] - 
                    mm[1:(N-lag),,drop=F])/mm[1:(N-lag),,drop=F]
   if (vec) pchange <- pchange[,1]
   tfclass(pchange) <- cls
   if (!is.null(tf)) tframed(pchange, tf) else pchange
}

percent.change.TSestModel <- function(model, base=NULL, lag=1, cumulate=F, e=F)
  {#The percent change calculation is done
   # with $estimates$pred and the result is an object of class TSdata
   TSdata(output=percent.change(model$estimates$pred))
  }

percent.change.TSdata <- function(data, base=NULL, lag=1, cumulate=F, e=F)
  {# The percent change calculation is done
   # with input and output and the result is an object of class TSdata.
   if (0 != (input.dimension(data)))  input.data(data)  <- percent.change(input.data(data))
   if (0 != (output.dimension(data))) output.data(data) <- percent.change(output.data(data))
   data
  }



# standardize <- function(ser){ # old version for a vector
#   m <- mean(ser)
#   v <- var(ser-m)
#   (ser-m)/(v^0.5)
#}

standardize <- function(ser){
	if (!is.matrix(ser)) stop("series should be a matrix.")
	means <- apply(ser, 2, mean)
	new <- ser - t(matrix(means, ncol(ser), nrow(ser)))
	svd.cov <- svd(var(new))
        scalefac <- svd.cov$u %*% diag(1/svd.cov$d^0.5, ncol = length(svd.cov$d))
        new <- new %*% t(scalefac)
	tframed(new, tf=tframe(ser), names=series.names(ser))
    }


############################################################

#     model and data scaling functions   <<<<<<<<<<

############################################################

# the following makes a generic copy of the function in the library to resolve
#  problems with version of S which do not have a generic function.
#  ( see also .First.lib at the beginning of dse1a.s
 # otherwise the following produces warning messages
invisible(
#   if (exists("scale.default"))
#     {scale.default <- scale.default
#      scale <- scale
#     }
   if (!exists("scale.default"))
     {if (exists("scale")) scale.default <- scale  
      #scale <- function(obj, ...) UseMethod("scale")
      scale <- function(x, ..., scale = TRUE) UseMethod("scale")
     }
  )
  


scale.TSdata <- function(data, scale) 
{# scale should be a list with two matrices or vectors, named input and output,
 # giving the multiplication factor for inputs and outputs.
 # vectors are treated as diagonal matrices.
 # If input or output are NULL then no transformation is applied.
 # The resulting data has inputs and outputs which are different from
 #  the original in proportion to scale. ie. if S and T are output and input
 #  scaling matrices then 
 #          y'(t) = S y(t) where y' is the new output
 #          u'(t) = S u(t) where u' is the new input
 sc <- input.data(scale)
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("input scale must be a vector or matrix")
    d <- input.data(data)
    tf <- tframe(d)
    names <- series.names(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    series.names(d) <- names
    input.data(data) <- d 
   }
 sc <- output.data(scale)
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("output scale must be a vector or matrix")
    d <- output.data(data)
    tf <- tframe(d)
    names <- series.names(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    series.names(d) <- names
    output.data(data) <- d 
   }
 data
}


scale.TSestModel <- function(model, scale) {scale(TSmodel(model), scale)}

scale.check <- function(obj, ...) UseMethod("scale.check")

scale.check.TSestModel <- function(model, scale){scale.check(TSmodel(model), scale)}

scale.check.TSmodel <- function(model, scale) 
{# This function only checks for some error conditions.
 if (!is.null(input.data(scale)))
   {if (is.matrix(input.data(scale)))
       {if (any(svd(input.data(scale))$d == 0))  
        stop("input.data(scale) transformations must be non singular.")
       }
    else if(any(input.data(scale)== 0)) stop("input.data(scale) elements must be non zero.")
   }
 if (!is.null(output.data(scale)))
   {if (is.matrix(output.data(scale)))
      {if (any(svd(output.data(scale))$d == 0))  
       stop("output.data(scale) transformations must be non singular.")
      }
    else if(any(output.data(scale)==0)) stop("output.data(scale) elements must be non zero.")
   }
 invisible(T)
}

# scale.SS <- function(model, scale){scale.check(model, scale)}

scale.innov <- function(model, scale)
{if (!scale.check(model, scale)) stop("scaling error.")
 if (!is.null(input.data(scale))) 
   {sc <- input.data(scale)
    if(is.vector(sc)) sc <- diag(sc, input.dimension(model))
    model$G <- model$G %*% solve(sc)
   }
 sc <- output.data(scale)
 if(is.vector(sc)) sc <- diag(sc, output.dimension(model))          
 model$H <- sc %*% model$H
 model$K <- model$K %*% solve(sc)
 set.parameters(model)
}

scale.non.innov <- function(model, scale)
{if (!scale.check(model, scale)) stop("scaling error.")
 if (!is.null(input.data(scale))) 
   {sc <- input.data(scale)
    if(is.vector(sc)) sc <- diag(sc, input.dimension(model))
    model$G <- model$G %*% solve(sc)
   }
 sc <- output.data(scale)
 if(is.vector(sc)) sc <- diag(sc, output.dimension(model))          
 model$H <- sc %*% model$H
 model$R <- sc %*% model$R %*% solve(sc)
 set.parameters(model)
}

scale.ARMA <- function(model, scale)
{if (!scale.check(model, scale)) stop("scaling error.")
 sc <- output.data(scale)
 if(is.vector(sc)) sc <- diag(sc, output.dimension(model))   
 model$A <- polyprod(sc, polyprod(model$A, solve(sc)))
 model$B <- polyprod(sc, polyprod(model$B, solve(sc)))
 if (!is.null(model$C)) 
   {sci <- input.data(scale)
    if (!is.null(sci)) 
       {if(is.vector(sci)) sci <- diag(sci, input.dimension(model))          
        model$C <- polyprod(model$C, solve(sci))
       }
    model$C <- polyprod(sc, model$C)
   }
 if (!is.null(model$TREND))  model$TREND <- sc %*% model$TREND
 set.parameters(model)
}


#######################################################################

#    test functions for dse1a.s dse1b.s dse1c.s and dse1d.s   <<<<<<<<<<

#######################################################################



dse1.function.tests <- function( verbose=T, synopsis=T, fuzz.small=1e-14, fuzz.large=1e-10)
{# A short set of tests of the main DSE functions using 
 #    eg1.DSE.data.diff.
 # The main short coming of these tests is that they do not test
 # functions which produce output, such as display, summary, graph
 # and check.residuals.
 # Note- using total sample.   Working paper estimated with sub-sample!

#  if (verbose) cat("dse1 test 7 ... ")
#  ok <- McMillan.degree(VARmodel$model, verbose=F)$distinct ==
#                McMillan.degree(ARMAmodel, verbose=F)$distinct
#  all.ok <- all.ok & ok 
#  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (!is.TSdata(eg1.DSE.data.diff))
     stop("Test data not found. Testing stopped.\n")
  max.error <- NA
  if (synopsis & !verbose) cat("All dse1 (kernel) tests ...")

  if (verbose) cat("dse1 test 0 ... ")
  # check "window"
  z <- tfwindow(output.data(eg1.DSE.data.diff), start=c(1980,1), end=c(1980,1))
  ok <- all( c (c(1,3)==dim(z), c(1980,1)==start(z), c(1980,1)==end(z)))
  z <- tfwindow(output.data(eg1.DSE.data.diff), start=c(1980,1), end=c(1982,12))
  ok <- ok & all( c (c(36,3)==dim(z), c(1980,1)==start(z), c(1982,12)==end(z)))
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 1 ... ")
  z <- est.VARX.ls(eg1.DSE.data.diff)
#  z <-eg1.DSE.data.diff
#  lsfit produces warning messages in the following
#  z$output[100,] <-NA
#  z <- est.VARX.ls(z, warn=F)
  VARmodel  <-  est.VARX.ar(eg1.DSE.data.diff, re.add.means=F, warn=F)
  SSmodel  <- to.SS(VARmodel)
  ok <- fuzz.large > abs(VARmodel$estimates$like[1] -
               l(SSmodel, eg1.DSE.data.diff, warn=F)$estimates$like[1])
  ok <- ok & is.TSestModel(VARmodel) & is.TSmodel(VARmodel$model)
  ok <- ok & (input.dimension(VARmodel) == input.dimension(SSmodel))
  ok <- ok & (input.dimension(VARmodel) == input.dimension(VARmodel$data))
  ok <- ok & (output.dimension(VARmodel) == output.dimension(SSmodel))
  ok <- ok & (output.dimension(VARmodel) == output.dimension(VARmodel$data))
  VARmodelB <- TSmodel(VARmodel)
  B <- t(chol(VARmodel$estimates$cov))
  VARmodelB$B <- array(B, c(1,dim(B)))  # has B != I
  VARmodelB <- set.parameters(VARmodelB)
  VARmodelB <- l(VARmodelB,VARmodel$data, warn=F)
  error <- max(abs(VARmodel$estimates$pred -VARmodelB$estimates$pred))
  ok <- ok & fuzz.large > error 
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 2 ... ")
  error <- abs(VARmodel$estimates$like[1] -
    l(set.arrays(SSmodel), eg1.DSE.data.diff,warn=F)$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 3 ... ")
  error <- abs(VARmodel$estimates$like[1]  - l(set.arrays(VARmodel), 
                          eg1.DSE.data.diff, warn=F)$estimates$like[1])
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 4 ... ")
  ARMAmodel <- to.ARMA(SSmodel)
  error <- abs(VARmodel$estimates$like[1] -
             l(ARMAmodel, eg1.DSE.data.diff, warn=F)$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 5 ... ")
  error <- abs(VARmodel$estimates$like[1] -
            l(ARMAmodel, eg1.DSE.data.diff,warn=F)$estimates$like[1])
  ok <- fuzz.large > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 6 ... ")
  error <- max(abs(sort(Mod(roots(TSmodel(VARmodel),by.poly=T))) -
                            sort(Mod(roots(SSmodel))) ))
  ok <-      fuzz.small > error
  err2 <- max(abs(sort(Mod(roots(TSmodel(VARmodel),by.poly=F))) -
                            sort(Mod(roots(SSmodel))) ))
  ok <- ok & fuzz.small > err2
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error,err2)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 7 ... ")
  d <-20
  true.roots <- c(-1/seq(d),1/seq(d),-seq(d),seq(d)) 
  A <- array(0, c(2,length(true.roots),length(true.roots)))
  A[1,,] <- diag(1,length(true.roots))
  A[2,,] <- diag(-true.roots, length(true.roots))
  # the following relies on roots using by.poly=F
  if(is.Splus()) options(expressions=1024)
  error <- max(Mod(
       sort(roots( ARMA(A=A, B=diag(1,length(true.roots)) ),by.poly=F))
     - sort(true.roots)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  #  complex roots test

  i <- pi*(1:10)/10.1  # this is half circle, but conjs also get generated
  # div by 10.1 instead of 10 prevents a real root with = conj
  true.roots <- complex(real=cos(i), imaginary=sin(i))  # on unit circle
  # scale simplifies sorting
  true.roots <- c(true.roots*(1+.2*1:10), true.roots*(1:10)/10) 
  A <- array(0, c(3,length(true.roots),length(true.roots)))
  A[1,,] <- diag(1,length(true.roots))
  A[2,,] <- diag(-2*Re(true.roots), length(true.roots))
  A[3,,] <- diag(Re(true.roots*Conj(true.roots)), length(true.roots))
  est.roots <- roots( ARMA(A=A, B=diag(1,length(true.roots)) ))
  ec <- 0<=Im(est.roots)
  error <- max(Mod(est.roots[ ec][order(Mod(est.roots[ ec]))]
                         - true.roots[order(Mod(true.roots))]))
  ok <- ok & fuzz.small > error
  err2 <- max(Mod(est.roots[!ec][order(Mod(est.roots[!ec]))]
                     - Conj(true.roots)[order(Mod(true.roots))]))
  ok <- ok & (fuzz.small > err2)
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error,err2)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 8 ... ")

  z  <- simulate(TSmodel(VARmodel), input=input.data(eg1.DSE.data.diff)) 
  zz <- simulate(TSmodel(VARmodel), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff))
  ok <- test.equal(z, zz, fuzz=fuzz.small)

  sigma <- solve(t(VARmodelB$model$B[1,,]) %*% VARmodelB$model$B[1,,])
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zzz <- simulate(TSmodel(VARmodelB), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff), SIGMA=sigma)
  error <- max(abs(output.data(z) - output.data(zzz)))
  ok <- ok & test.equal(z, zzz, fuzz=fuzz.small)

  # next use estimates$cov
  z  <- simulate(VARmodel, input=input.data(eg1.DSE.data.diff)) 
  sigma <- VARmodel$estimates$cov
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zz <- simulate(TSmodel(VARmodel), rng=get.RNG(z), 
                     input=input.data(eg1.DSE.data.diff), SIGMA=sigma)
  ok <- ok & test.equal(z, zz, fuzz=fuzz.small)

  ok <- ok & test.equal(summary(zz)$estimates,
                        summary(zz)$estimates, fuzz=fuzz.small)

  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse1 test 9 ... ")
  z  <- simulate(SSmodel, input=input.data(eg1.DSE.data.diff)) 
  ok <- test.equal(z,simulate(SSmodel, rng=get.RNG(z), 
                      input=input.data(eg1.DSE.data.diff)))
  ok <- ok & test.equal(summary(z)$estimates,
                        summary(z)$estimates, fuzz=fuzz.small)

  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse1 test 10... ")
  ok <- stability(SSmodel, verbose=F)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse1 test 11 ...")
  scale.fac <- diag(1:3)
  scale.fac[1,3] <-.5
  scale.pred <- VARmodel$estimates$pred %*% t(scale.fac)
  scale.fac <- list(output=scale.fac)
  error <- max(abs(scale.pred -
        l(scale(VARmodel$model, scale.fac), 
          scale(eg1.DSE.data.diff, scale.fac), warn=F)$estimates$pred))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
  if (!ok)
    cat("scale in dse1 test 11 causes an error if DSE is not at the beginning of the search() list")

  if (verbose) cat("dse1 test 12... ")
  error <- max(abs(scale.pred
         - l(scale(SSmodel, scale.fac), 
             scale(eg1.DSE.data.diff, scale.fac))$estimates$pred))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse1 test 13... ")
  z <- eg1.DSE.data.diff
  ok <- test.equal(z,
      TSdata(output=output.data(combine(z,z), series=seq(output.dimension(z))),
              input= input.data(combine(z,z), series=seq( input.dimension(z))))) 
 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All dse1 (kernel) tests completed")
     if (all.ok) cat(" OK\n")
     else cat(", some FAILED! max.error = ", max.error,"\n")
    }
  invisible(all.ok)
}


#######################################################################

#                    end

#######################################################################


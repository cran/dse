#   2000/04/18 11:15:50 
#######################################################################

# These functions use functions in curvature.hs

#######################################################################

#                            curvature calculation

#######################################################################

curvature.TSestModel <- function (emodel, warn=T, ...)  
 {curvature(genD(emodel, ...),  warn=warn)}

#-----------------------------------------------------------------------

# S routines for calulating curvatures a la Bates and Watts.

# Notes:
#   Generating the D matrix can be computationally demanding. There are three
#   (or four) versions of the code for generating this matrix. The S version is slow 
#   and suffers from memory problems due to a bug in the way the current version of S
#   (S-PLUS Version 3.0 for SUN4) allocates memory in loops. The C version is
#   faster but suffers (even worse) from the memory problem. (A fix is promised 
#   in the next release of S-plus.) Both the S and the C versions take the name
#   of an S function as an arguement and call it. The compiled version is fast and
#   does not suffer from the memory problem, but works only with ARMA 
#   and KF models.
#    
#-----------------------------------------------------------------------

#######################################################################

#                               D matrix calculation

#######################################################################

genD.TSestModel <- function(estModel, ...)
   { invisible(genD( TSmodel(estModel), TSdata(estModel), ...))}

#genD.TSmodel <- function(model, data, ...) {NextMethod("genD.TSmodel")}

genD.ARMA <- function(model, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   n <-length(c(output.data(data))) # this has the same length as the residual
   sampleT <-periods(data) 
   ps <-dim(model$A)[3]
      ms <-dim(model$C)[3]
      ns <-ps  # this could be 1 except z0 is used for TREND and F, Q and R for scratch
     loc   <- match(model$location,       c("A","B","C","t"))
     cloc  <- match(model$const.location, c("A","B","C","t"))
   if(is.null(ms))
     {ms<-0
      C<-array(0,c(1,ns,1)) # can't call compiled with 0 length arrays
      u <- matrix(0,sampleT,1)
      G <- matrix(0,ns,1)
     }
   else
     {C <-model$C
      u <- input.data(data)
      G <- matrix(0,ns,ms)
     } 
   x <-model$parms
   zt <- 0.0001
   #Rbug zt <-unix.time(l(model,data,sampleT,sampleT,result="like"))
   zt <- zt[1] *r*2*(length(x)*(length(x) + 3))/2
   if(zt>30) 
       cat("D matrix calculation very roughly estimated to take about ",
            (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(1e20, n,(length(x)*(length(x) + 3))/2) 
   Daprox <-  matrix(0,n,r) 
   Hdiag  <-  matrix(0,n,length(x))
   Haprox <-  matrix(0,n,r)
   storage.mode(D)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   D <-.Fortran("gend",
            D=D,
            as.integer(is.ARMA(model)), 
            p=as.integer(length(x)),
            x=as.double(x),
            as.double(h0),
            as.integer(n),    #6
            as.integer((length(x)*(length(x)+ 3))/2), #cols of D
            f0=double(n),      
            r=as.integer(r),
            #                       work space for GEND
            Haprox=Haprox,     #10   
            Hdiag=Hdiag,         
            Daprox=Daprox,        
            x=double(length(x)),
            delta=double(length(x)),
            f1=double(n),
            f2=double(n),
            #                   work space for ARMAp / KFp
         #     cov=matrix(1e20,ps,ps),       
            # pred is f0,f1,f2 passed above
            as.integer(ms),     #  input dimension m  #17
            as.integer(ps),     # output dimension p 
            as.integer(sampleT),   
            as.integer(periods(data)), 
            u=as.double(u), 
            y=as.double(output.data(data)),   
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location), #23
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=as.double(model$const),
            as.integer(cloc),  #as.character(model$const.location),
            as.integer(model$const.i ), 
            as.integer(model$const.j),
        #  for ARMA models:
            as.integer(model$l),       #31
            as.integer(model$const.l),
            as.integer( dim(model$A)[1]),#1+order of A
            as.integer( dim(model$B)[1]),#1+order of B
            as.integer( dim(C)[1]),#1+order of C
            A=model$A,   
            B=model$B,  
            C=C, 
        #  for state space models:
            as.integer(ns),  # state dim.  #39
        #    state=as.double(matrix(0,sampleT,ns)),  
        #    track=as.double(array(0,c(sampleT,ns,ns))),  
            z0=as.double(rep(0,ps)), # note: this is TREND for ARMA
            P0=as.double(diag(0,ps)),
            F=as.double(matrix(0,ns,ns)),  
            G=as.double(G),  
            H=as.double(matrix(0,ps,ns)),  
            K=as.double(matrix(0,ns,ps)),  
            Q=as.double(matrix(0,ns,ns)),  
            R=as.double(matrix(0,ps,ps)), 
            gain=as.logical(F)   #48
            )[c("D","p","f0", "x", "r")]
   D$d   <- d
   D$eps <- eps
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(output.data(data))
   invisible(classed(D,"Darray")) #constructor
}

#Rbug this does not seem to work as genD.KF or genD.KF.innov

genD.innov <- function(model, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   n <-length(c(output.data(data))) # this has the same length as the residual
   sampleT <-periods(data) 
   ns <-dim(model$F)[2]
   ps <-dim(model$H)[1]
   ms <-dim(model$G)[2]
   if(is.null(ms))
     {ms<-0
      C<-array(0,c(1,ns,1)) # can't call compiled with 0 length arrays
      u <- matrix(0,sampleT,1)
      G <- matrix(0,ns,1)
     }
   else
     {C <- array(0,c(1,ns,ms))
      u <- input.data(data)
      G <- matrix(0,ns,ms)
     } 
   x <-model$parms
   zt <- 0.0001
   #Rbug zt <-unix.time(l(model,data,sampleT,sampleT,result="like"))
   zt <- zt[1] *r*2*(length(x)*(length(x) + 3))/2
   if(zt>30) 
       cat("D matrix calculation very roughly estimated to take about ",
           (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(1e20, n,(length(x)*(length(x) + 3))/2) 
   Daprox <-  matrix(0,n,r) 
   Hdiag  <-  matrix(0,n,length(x))
   Haprox <-  matrix(0,n,r)
   storage.mode(D)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   loc   <- match(model$location,       c("f","G","H","K","Q","R","z","P"))
   cloc  <- match(model$const.location, c("f","G","H","K","Q","R","z","P"))
   D <-.Fortran("gend",
            D=D,
            as.integer(is.ARMA(model)), 
            p=as.integer(length(x)),
            x=as.double(x),
            as.double(h0),
            as.integer(n),    #6
            as.integer((length(x)*(length(x)+ 3))/2), #cols of D
            f0=double(n),      
            r=as.integer(r),
            #                       work space for GEND
            Haprox=Haprox,        
            Hdiag=Hdiag,         
            Daprox=Daprox,        
            x=double(length(x)),
            delta=double(length(x)),
            f1=double(n),                # 15
            f2=double(n),
            #                   work space for ARMAp / KFp       
            # cov=matrix(1e20,ps,ps),       
            # pred is f0,f1,f2 passed above
            as.integer(ms),     #  input dimension m  
            as.integer(ps),     # output dimension p 
            as.integer(sampleT),   
            as.integer(periods(data)), 
            u=as.double(u), 
            y=as.double(output.data(data)),   #22
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location) bug
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=as.double(model$const),     #28
            as.integer(cloc),  #as.character(model$const.location),
            as.integer(model$const.i ), 
            as.integer(model$const.j),
        #  for ARMA models:
            as.integer(loc),                   #31
            as.integer(loc),
            as.integer(1),#1+order of A
            as.integer(1),#1+order of B
            as.integer(1),#1+order of C
            A=as.double(array(0,c(1,ps,ps))),   
            B=as.double(array(0,c(1,ps,ps))),  
            C=as.double(C), 
        #  for state space models:
            as.integer(ns),  # state dim.     #39
        #    state=as.double(matrix(0,sampleT,ns)),  
        #    track=as.double(array(0,c(sampleT,ns,ns))),  
            z0=as.double(rep(0,ns)),  
            P0=as.double(diag(0,ps)),
            F=as.double(matrix(0,ns,ns)),  
            G=as.double(G),  
            H=as.double(matrix(0,ps,ns)),  
            K=as.double(matrix(0,ns,ps)),  
            Q=as.double(matrix(0,ns,ns)),  
            R=as.double(matrix(0,ps,ps)),
            gain=as.logical(is.innov.SS(model)) #48
            )[c("D","p","f0", "x", "r")]
   D$d   <- d
   D$eps <- eps
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(output.data(data))
   invisible(classed(D, "Darray")) #constructor
}

#######################################################################

#              span (dimension of tangent space) calculation

#######################################################################

span.TSestModel <- function (emodel, compiled=.DSECOMPILED, ...)  
 {# calculate the singular values of the tangents
  # the compiled version calculates the whole D matrix (which seems like
  # a waste, but the compiled code is much faster, so ...
  if (compiled)
   {D <- genD(emodel, ... )$D[,seq(length(parms(emodel))),drop=F]
    return(svd(D)$d)
   }
  else {
     # there should be a better way to do this
     global.assign("func.residual.TSestModel",  function(parms,Shape,data)
      {c(l(set.arrays(Shape,parms=parms),data,result="pred") - output.data(data))})
     on.exit(remove("func.residual.TSestModel"))
     span.default("func.residual.TSestModel", parms(emodel),
               func.args=list(Shape=TSmodel(emodel), data=TSdata(emodel), ...))
     }
 }


#######################################################################

#            calculate Fisher Info (Hessian of the likelihood)

#######################################################################

hessian.TSestModel <- function (emodel, ...)  
 {# like returns neg. log likelihood
  # there should be a better way to do this
  global.assign("func.hessian.TSestModel",  function(parms,Shape,data)
  	{l(set.arrays(Shape,parms=parms),data,result="like")})
  on.exit(rm(func.hessian.TSestModel))
  hessian.default( "func.hessian.TSestModel", parms(emodel),
   func.args=list(Shape=TSmodel(emodel), data=TSdata(emodel), ...))
 }

# was func.args=append(list(Shape=TSmodel(emodel), data=TSdata(emodel)), list(...)))
#######################################################################

# Test routines for calulating curvatures.

#######################################################################
 


dsecurvature.function.testsA <- function( verbose=T, synopsis=T, fuzz.small=1e-12, 
    fuzz.large=1e-6, fuzz.very.large=1e-2, show.details=F, show.extra.details=F,
    print.values=T, digits=18)
{# Tests of DSE curvature functions

# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to RNG when R 1.0.0 was released
#   (Feb. 29, 2000) and also previously.

  random.number.test()
  
 test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

 all.ok <- list(ok=T, error=NA)

test <- function(tst, good, all.ok, flag="", fuzz=1e-14, print.values=F)
  {# initialize all.ok <- list(ok=T, error=NA)
   error <- max(abs(good-tst))
   if (any(is.na(error)) || any(is.nan(error)) || fuzz < error) 
     {if (print.values) print.test.value(c(tst), digits=18)
      all.ok$error <- if(is.na(all.ok$error)) error else max(error, all.ok$error)
      all.ok$ok <- F  # assumes ok is initialized T for this to work 
      cat(" ",flag, "failed! error = ", error,"\n")
     }
   all.ok
  }

  if (synopsis & !verbose) cat("All curvature tests A ...")
  if (verbose) cat("DSE curvature test A 1 ...")

# simplified from user guide
  VARmodel<-ARMA(A=array(c(1,.5,.3),c(3,1,1)),
                 B=array(1,c(1,1,1)),
                 C=NULL, description="simplified guide example")
  VARmodel <-l(VARmodel,simulate(VARmodel, rng=test.rng))

# unstable model
  VARmodel2<-ARMA(A=array(c(1,-0.5,-0.5),c(3,1,1)),
                 B=array(1,c(1,1,1)),
                 C=NULL, description="simplified guide example")
  VARmodel2 <-l(VARmodel2,simulate(VARmodel2, rng=test.rng))
# Mod(roots(VARmodel, by.poly=T))
# [1] 0.5477226 0.5477226
# Mod(roots(VARmodel2, by.poly=T))
# [1] 1.0 0.5

  if (fuzz.small < abs(VARmodel$estimates$like[1] - 143.78939695547763 ))
    {warning("VARmodel  likelihood  does not correspond to expected value.")
     cat("VARmodel  likelihood:")
     print(VARmodel$estimates$like[1], digits=digits)
    }
 
  if( fuzz.small < abs(sum(VARmodel$data$noise$w) - 9.384686064093962  ))
    {warning("Check sum of noise does not correspond to expected test value.")
     cat("Check sum of noise:")
     print(sum(VARmodel$data$noise$w), digits=digits)
    }

  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)
  if (fuzz.small < abs(ARMAmodel$estimates$like[1]- VARmodel$estimates$like[1]))
    {warning("ARMAmodel likelihood does not correspond to expected test value.")
     cat("ARMAmodel  likelihood:")
     print(ARMAmodel$estimates$like[1], digits=digits)
    }

  spanVAR <- span(VARmodel, compiled=F)

  all.ok <- test(spanVAR, c(12.583380392358416,  9.21209438442244277),
                 all.ok, flag="1a", fuzz=10*fuzz.small, #10* for Linux vs Solaris R
                 print.values=print.values)
 
  spanVAR.f <- span(VARmodel, compiled=.DSECOMPILED)

  all.ok <- test(spanVAR.f, spanVAR, all.ok, flag="1b",  fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 2 ...")

  spanSS.f <- span(SSmodel, compiled=.DSECOMPILED)

  tst <- c( 11.5173146571970904,   9.8419970117338984,
             5.42943770603545595,  1.65151435950942771)

  all.ok <- test(spanSS.f, tst, all.ok, flag="2a", fuzz=fuzz.large, #fuzz.small works in Solaris
                 print.values=print.values)
  spanSS <- span(SSmodel, compiled=F)
  all.ok <- test(spanSS, spanSS.f, all.ok, flag="2b", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 3 ...")
  spanARMA.f <- span(ARMAmodel, compiled=.DSECOMPILED)

  all.ok <- test(spanARMA.f, spanVAR, all.ok, flag="3a", fuzz=fuzz.small,
                 print.values=print.values)


  spanARMA <- span(ARMAmodel, compiled=F)
  all.ok <- test(spanARMA, spanARMA.f, all.ok, flag="3b", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 4 ...")
# calculating from the score using genD might be quicker.
  hessianVAR <- hessian(VARmodel) 

  tst <- c( 9.83444900505639374,  -5.64196107342148512,
           -5.64196107342148512,  116.345053001385807)

  all.ok <- test(hessianVAR, tst, all.ok, flag="4a", fuzz=fuzz.large,
                 print.values=print.values)


  hessianSS <- hessian(SSmodel)

  tst <-  c( 3.95368029144797717,  -0.418007621718194833,  -9.41644138321294832,
            36.9034481180689653,  -0.418007621718194833,  1.68828078196901132,  
	    39.9657123856431085,  -16.2517115170908539,  -9.41644138321294832,  
	    39.9657123856431085,  10.2620303429070088,  -65.7865038150653731,  
	    36.9034481180689653,  -16.2517115170908539,  -65.7865038150653731,  
	   107.778335544667868)

  all.ok <- test(hessianSS, tst, all.ok, flag="4b", fuzz=10*fuzz.large,
                 print.values=print.values)


  hessianARMA <- hessian(ARMAmodel)
  all.ok <- test(hessianARMA, hessianVAR, all.ok, flag="4c", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 5 ...")
  curvatureVAR <- curvature(VARmodel)$stats
  all.ok <- test(curvatureVAR, c(2, 100, 0.05, 0, 0, 0, 0, 1, 1), all.ok, flag="5a", fuzz=fuzz.large,
                 print.values=print.values)

# and for comparison
# there should be a better way to do this
  global.assign("func.residual",  function(parms,Shape,data)
    {c(l(set.arrays(Shape,parms=parms),data,result="pred") - output.data(data))})
  on.exit(rm(func.residual))
     
  curvatureVAR.def <- curvature.default("func.residual", parms(VARmodel), 
                  func.args=list(Shape=TSmodel(VARmodel),data=TSdata(VARmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=F)$stats
  all.ok <- test(curvatureVAR.def, curvatureVAR, all.ok, flag="5b", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 6 ...")
  curvatureSS <- curvature(SSmodel, warn=F)$stats
  # neg sqrt in axis ratio produces warning if warn=T

# previous rng
#  tst <- c(4, 100, 0.05,  0.9992455609137233,  0.9470084226509906, 
#     1.56931709092204,  1.487278563994181,  1.000000000828572, 
#     1.321413857107741)
 
  tst <- c(4, 100, 0.05, 1.02548364377805989,  0.979897817007443384,  
     1.61052405093533291,  1.53893142160253493,  1.00000000114714904)  #, NaN

  all.ok <- test(curvatureSS[-9], tst, all.ok, flag="6a", fuzz=fuzz.large,
                 print.values=print.values)


  curvatureSS.def <- curvature.default("func.residual", parms(SSmodel), 
               func.args=list(Shape=TSmodel(SSmodel), data=TSdata(SSmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=F, warn=F)$stats
  # neg sqrt in axis ratio produces warning if warn=T

  all.ok <- test(curvatureSS.def[-9], curvatureSS[-9], all.ok, flag="6b", fuzz=10*fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 7 ...")
  curvatureARMA <- curvature(ARMAmodel)$stats
  all.ok <- test(curvatureARMA, curvatureVAR, all.ok, flag="7a", fuzz=fuzz.small,
                 print.values=print.values)
 
  curvatureARMA.def <- curvature.default("func.residual", parms(ARMAmodel), 
               func.args=list(Shape=TSmodel(ARMAmodel), data=TSdata(ARMAmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=F)$stats

  all.ok <- test(curvatureARMA.def, curvatureARMA, all.ok, flag="7b", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 8 ...")
# from user guide

  VARmodel<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 

# Note this gives a terrible fit.
  VARmodel<-l(VARmodel,simulate(VARmodel, rng=test.rng))
  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)

  spanVAR.f <- span(VARmodel, compiled=.DSECOMPILED) 

#  if(is.Splus()) tst <-  c(
#      18.28342261104111799, 16.86252848720547703, 15.47359559657559025,
#      12.60857465926911836, 10.85404071073360299, 10.80245325753396202,
#       9.79091689469922066,  9.21160820465702379,  2.32881291213965458,
#       1.97834187456894184,  1.50902538655984464 )
#  if(is.R())     tst <-  c(
#        17.8366876297672725, 16.7671374021287249, 16.2282896291638892,
#        13.3626055311535605, 10.3267790597707432,  9.9307009769349008,
#         8.8813776725358053,  8.2133988282015906,  2.2202966179396091,
#         1.9141400784167535,  1.4838767428262707 )     #f77 on Sun5 R0.61.1
#     c( 1.783668762976735e+01, 1.676713740212845e+01, 1.622828962916468e+01,
#     1.336260553115620e+01, 1.032677905977513e+01, 9.930700976937008,
#     8.881377672535175, 8.213398828201516, 2.220296617939313,
#     1.914140078418723, 1.483876742824580)))  #R 0.49
#  tst <-c(13.30038806708721,  12.23076792331929,  9.752849398484228, 
#         8.786241709293281,  8.583699014891037,  7.624249106943636,  
#         5.998914711634687,   5.61060655287252,  1.774590614465879,  
#         1.436530005339643,  0.9712521284803599)
   tst <- c( 17.3672899077923439,  15.9424331056250885,  14.5699882827285894,
             11.6462671196918794,  10.7723733430260129,  10.5668728474602247, 
	      8.76191214115850769,  7.93300083097722553,  2.24655970809629668,  
	      1.88819440545349604,  1.47455500924911109)
  all.ok <- test(spanVAR.f, tst, all.ok, flag="8a", fuzz=fuzz.large,
                 print.values=print.values)
 
  spanVAR <- span(VARmodel, compiled=F) 
  all.ok <- test(spanVAR, spanVAR.f, all.ok, flag="8b", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 9 ...")
  spanSS <- span(SSmodel, compiled=.DSECOMPILED)

#  if(is.Splus()) tst <-  c(
#         12.98465581800713764, 12.04846041797238776, 11.75207879184403659,
#         10.86851614803762800,  9.41015226452319986,  9.14102587417263912,
#         8.81980414072582874,  8.32921589456761780,   5.12046284222839354,
#         5.02836144306142163,  3.44776911209970471,   3.14630978631044567,
#         2.20492937840385128,  2.07975512652550254,   1.93518614772789666,
#         1.87713878204627571 )
#  if(is.R())     tst <-  c(
#         14.3191707820078413, 12.7553932264152472, 12.0399559044862077,
#         11.1780750248459100,  8.7642374743163280,  8.4110650474333628,
#          8.1119238141018251,  7.7112054454257013,  4.5816849188761744,
#          4.4928124488652665,  3.0909505682400633,  2.7790878709102365,
#          2.2901166620448916,  2.1211736442161429,  1.8231113211019712,
#          1.7364603432999757 )         #f77 on Sun5 R0.61.1

#      c( 1.431917078200739e+01, 1.275539322641331e+01, 1.203995590448420e+01,
#       1.117807502484729e+01, 8.764237474316120, 8.411065047433565,
#       8.111923814106442,     7.711205445427577, 4.581684918879841,
#       4.492812448867578,     3.090950568240565, 2.779087870909759,
#       2.290116662046349,     2.121173644211920, 1.823111321101024,
#       1.736460343300366)))  # 0.49

#  tst <-  c(9.470240142187398,  8.880447829525087,  7.552501916476404,
#            7.531576128069482,  7.169071366088333,  6.868663706072056,  
#            5.629811912135156,  5.343700909140345,  3.984898912643869,  
#            3.899294396457272,  2.444843033246666,  2.207333797038643,  
#            1.447912667824839,  1.363598776819928,  1.256253315103265,  
#            1.205277678709112)
  tst <- c(12.6684354906294967,  11.4927966118834544,  11.0537104821112262,
           10.3371110331134979,  9.77605886766110288,  9.21103083988342064, 
	    8.16488126052858476,   7.7872986620231881,  4.79640068843280787, 
	    4.69148925882611323,  3.23108014824838419,  2.94195878161703162,  
	    2.11511171224067507,  2.01921578796399359,  1.85314036275455907,  
	    1.79333516556799233)
  all.ok <- test(spanSS, tst, all.ok, flag="9a", fuzz=fuzz.large,
                 print.values=print.values)
  spanSS <- span(SSmodel, compiled=F) 
  all.ok <- test(spanSS, tst, all.ok, flag="9b", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 10...")
  spanARMA.f <- span(ARMAmodel, compiled=.DSECOMPILED)

# if(is.Splus()) tst <-  c(
#     2.77119952678264987e+01, 2.53919116411379555e+01, 2.32812418066877100e+01,
#     1.90837600070161528e+01, 1.73322831223985432e+01, 1.48319783197278330e+01,
#     1.41417499505086042e+01, 1.19871189922967893e+01, 1.16051570137965676e+01,
#     1.08623372662499182e+01, 7.21807567275184692e+00, 6.80015402524795309e+00,
#     5.99427285923686615e+00, 5.84925830647237177e+00, 3.02327766004096032e+00,
#     2.67973563219726740e+00, 7.50824755195420801e-01, 4.76770151276800724e-01,
#     6.30824258854855713e-15, 1.22839566408688817e-15, 6.66516531623418467e-17,
#      0.00000000000000000e+00 )
#  if(is.R())     tst <-  c(
#       3.0111693007925084e+01, 2.6903738538557757e+01, 2.1949166927772083e+01,
#       1.9412118834825545e+01, 1.7023849849248617e+01, 1.5059167312762057e+01,
#       1.2198080331280515e+01, 1.1203212614180140e+01, 1.0534424040150872e+01,
#       9.7611634784175791e+00, 6.7724158272123773e+00, 6.0772753814400602e+00,
#       5.5378918184900012e+00, 5.2255556591227439e+00, 2.7020830182124316e+00,
#       2.5433734102944965e+00, 7.4778345526592260e-01, 4.6199586045921071e-01,
#       5.5366294436852289e-07, 2.5570412409273409e-07, 7.3166444674331792e-17,
#       0.0000000000000000e+00)         #f77 on Sun5 R0.61.1

#   c( 3.011169300794247e+01, 2.690373853859740e+01, 2.194916692782002e+01,
#      1.941211883487011e+01, 1.702384984930141e+01, 1.505916731276925e+01,
#      1.219808033127580e+01, 1.120321261421515e+01, 1.053442404016652e+01,
#      9.761163478425873e+00, 6.772415827211282e+00, 6.077275381439208e+00,
#      5.537891818483530e+00, 5.225555659106909e+00, 2.702083018243574e+00,
#      2.543373410280471e+00, 7.477834552669831e-01, 4.619958604596496e-01,
#      2.038270252146126e-06, 2.556984253017427e-07, 1.795257106881694e-15,
#      7.973651962079763e-17))) #R0.49
#      c( 2.768304162846298e+01, 2.393941491091734e+01, 1.988022205000436e+01, 
#         1.822550017859639e+01, 1.630502757392022e+01, 1.461512246650834e+01, 
#         1.275451631155268e+01, 1.221047501548486e+01, 1.073726354564602e+01, 
#         1.026204415323931e+01, 7.265742818227836e+00, 6.282013694588179e+00, 
#         5.974725282653903e+00, 5.251293695221127e+00, 2.860581954100438e+00, 
#         2.471994703685991e+00, 7.121798982572878e-01, 4.406741938800671e-01, 
#         1.433057563479431e-15, 9.026520997009618e-16, 4.800022955996391e-17, 
#         0.000000000000000e+00 )))

# R Solaris
  tst <-  c(3011.318489057963,  19.37947620116293,  17.47202335066703,
            15.20951506186483,  13.78393160485094,  12.03840540489752, 
             9.92428984973578,   9.697246859100282,  9.244229358342194,
             7.152271025275832,  7.092732858677983,  5.69297886809311,  
             4.736206089454028,  4.304104324296326,  3.708282982105869, 
             2.11735492156378,   1.95253799784708,   0.5057146032686363, 
             0.3160112544521056, 2.867793448120931e-15, 8.538596555411535e-16,
                  0)

# R Linux
  tst <- c( 19.5069136184747691,   17.473349776355132,  15.3942712009364495,
            13.7862136356176972,  12.0470502264922672,  9.93058117606233282,
            9.70738285213095686,  9.39682286923732946,  7.15244420534333436,
            7.14045582919181054,  5.72442802958446784,  4.74384059078413145,
            4.3178723214231276,   3.7093267481375447,  2.13283668011786354,
            1.95407955864386418,  0.506460713920498851,  0.316713386010944509,
            8.73275433068441726e-16,  2.91045244991302892e-16)

# Splus Solaris values:
#c(  [1]  2.854799065168644e+03  1.927325423818564e+01  1.745392647304185e+01
# [4]  1.508940688448259e+01  1.377654544670931e+01  1.203796836342551e+01
# [7]  9.926711794026920e+00  9.705558069753085e+00  9.214674640019336e+00
#[10]  7.151912789329529e+00  7.110891344071947e+00  5.676460671597718e+00
#[13]  4.724769204414272e+00  4.311928305576296e+00  3.708361783746251e+00
#[16]  2.130165310609959e+00  1.953704281613287e+00  5.057391241655013e-01
#[19]  3.137209572167103e-01  1.310169543169194e-15  5.690213922485003e-16
#[22]  0.000000000000000e+00 )
warning("Skipping 10a comparison. Problem is too ill-conditioned.")

#  all.ok <- test(spanARMA.f, tst, all.ok, flag="10a", fuzz=fuzz.large,
#                 print.values=print.values)
 
  spanARMA <- span(ARMAmodel, compiled=F) 
  all.ok <- test(spanARMA, spanARMA.f, all.ok, flag="10b", fuzz=fuzz.small,
                 print.values=print.values)
 

  ARMAmodel.fixed <- l(fix.constants(ARMAmodel), VARmodel$data)
  spanARMA.fix <- span(ARMAmodel.fixed)

#  if(is.Splus()) tst <-  c(
#     27.711995267845011881, 25.391911641137962619, 23.281241806692566598,
#     19.083760007020678984, 17.332283122405144127, 14.831978319714982817,
#     14.141749950539304947, 11.987118992298254838, 11.605157013791078668,
#     10.862337266239681099,  7.218075672755373873,  6.800154025249395495,
#      5.994272859236035700,  5.849258306460957790,  3.023277660033443670,
#      2.679735632202348672,  0.750824755198911564,  0.476770151275852538 )
#  if(is.R())     tst <-  c(
#     30.11169300792521142, 26.90373853855186326, 21.94916692777263734,
#     19.41211883482275269, 17.02384984926177935, 15.05916731276262510,
#     12.19808033128727054, 11.20321261418225234, 10.53442404015623701,
#      9.76116347841036180,  6.77241582721368829,  6.07727538144364043,
#      5.53789181847943901,  5.22555565912482223,  2.70208301821842412,
#      2.54337341031345465,  0.74778345526932322,  0.46199586045587798 ) 

#  tst <- c(19.50691361845315,  17.47334977636005,  15.39427120095003,
#          13.78621363562109,  12.04705022651371,  9.930581176044308,  
#          9.707382852080064,  9.396822869247679,  7.152444205356154,  
#          7.140455829182858,  5.724428029579137,  4.743840590786617,  
#          4.317872321424001,  3.709326748136878,  2.132836680138603,  
#          1.954079558657604,  0.5064607139184161,  0.3167133860055982)

  tst <- c(26.3659986552698982,  23.8359546246227616,  19.5634874084708699,
           17.2143457352275391,  15.1394364958679244,  14.7208570925179494,  
	   12.4221689849793311,  11.2407321096646733,  10.7899650394872069,  
	   10.1698037554743284,  7.99664809385088837,  6.32142067466109392,  
	    6.29099304634241019,  5.47437857637598668,  2.66712150086359934,  
	    2.51249349659001497,  0.683166443549642177,  0.459517858793767109)

  all.ok <- test(spanARMA.fix, tst, all.ok, flag="10c", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 11...")
  #  following test values have all been set using 
  #    R0.63.3pre and gnu f77 on SunOS 5.6 (Solaris)
  #  and reset as of R 1.0.0 
  hess <- hessian(VARmodel)  

# 1472.1174
  all.ok <- test(sum(hess), 1397.19588043396584, all.ok, flag="11a", fuzz=fuzz.very.large,
                 print.values=print.values)
 
  if (verbose) cat("b")
  hess <- hessian(SSmodel)
  #1393.953399146646
  all.ok <- test(sum(hess), 1351.52741934382425, all.ok, flag="11b", fuzz=fuzz.very.large,
                 print.values=print.values)
 
  if (verbose) cat("c")
  hess <- hessian(ARMAmodel) 
  # R Linux      -1409.96904582771185
  # R pre 1.0 Solaris   385220.9412876417
  # R 1.0.0   Solaris     -352.3883095147531
  # Splus          728.89769711477993
  #all.ok <- test(sum(diag(hess)), 385220.9412876417, all.ok, flag="11c", fuzz=fuzz.very.large,
  #               print.values=print.values)
  warning("Skipping 11c comparison. Problem is too ill-conditioned.")
  # print(sum(diag(hess)), digits=18)
 
  if (verbose) cat("d")
  hess <- hessian(ARMAmodel.fixed)
  # 3038.56
  all.ok <- test(sum(hess), 3039.23996170283772, all.ok, flag="11d", fuzz=fuzz.very.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 12...")
  curvatureVAR <- curvature(VARmodel,warn=F)$stats
#  tst <- c(11, 200, 0.05, 0.8490449316698463, 0.712275843318316,
#        1.151573691713139,  0.9660715137831059, 1.000000000576390) #, NaN)
  tst <- c(11, 200, 0.05, 0.807400747659583362,  0.681458988757835171,
         1.09509099576822266,  0.924274103952975268,  1.00000000053258109)
  all.ok <- test(curvatureVAR[-9], tst, all.ok, flag="12a", fuzz=fuzz.large,
                 print.values=print.values)
 
  curvatureVAR.def <- curvature.default("func.residual", parms(VARmodel), 
              func.args=list(Shape=TSmodel(VARmodel), data=TSdata(VARmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=F)$stats
  all.ok <- test(curvatureVAR.def[-9], curvatureVAR[-9], all.ok, flag="12b", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 13...")
  curvatureSS <- curvature(SSmodel, warn=F)$stats  
  # neg sqrt in axis ratio produces warning if warn=T
#  tst <- c(16, 200, 0.05, 1.492400232033412, 1.197222853812868,
#       1.945129668075371,   1.560408288784785,   1.000000000054351)
  tst <- c(16, 200, 0.05, 1.42377065952031878,  1.11780582522466698,
        1.85568086289751188,  1.45689958170905332,  1.00000000346571949)
  all.ok <- test(curvatureSS[-9], tst, all.ok, flag="13", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 14...")
  curvatureARMA.fixed <- curvature(ARMAmodel.fixed, warn=F)$stats
  # neg sqrt in axis ratio produces warning if warn=T
# R 0.64.1
#  tst <- c(18, 200, 0.05, 2.381219927481035, 2.270893359206098,
#            3.06872275232684,  2.92654283591333,  1.001605964896068)
# Splus 3.3
# tst <- c(18, 200, 0.05, 2.381219668018957, 2.270892355686405,
#           3.068722417953193, 2.926541542658688, 1.00160610746321)
# R 1.1.0 (devel)
  tst <- c(18, 200, 0.05, 2.34112169445768892,  2.29032831656791913,
            3.0170474078587568,  2.95158902973962212,  1.00047805155781822)
  all.ok <- test(curvatureARMA.fixed[-9], tst, all.ok, flag="14", fuzz=10*fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test A 15...")
 # next line max prob  and below too???
  curvatureARMA <- curvature(ARMAmodel, warn=F)$stats #singular matrix in solve(R11)

# D generated. Calculating curvature.Warning: eliminating degenerate subspace for R11.
# Warning: B is not symmetric.
# Warning: eigen for non-symmetric matrix is fictitious
#   With warn=T the following is generated:
# N.B. (I-B) is not positive definite as it should be at a local min! 
#    ARMAX will calculate to NA 
# Eigenvalues of B [1] -1.098230e+47 -5.735942e+15 -1.898362e+13 -5.863185e+01 -1.398720e+00
#  [6] -8.339023e-01 -2.658607e-01 -7.732492e-02 -6.136313e-02 -2.463111e-02
# [11] -2.361107e-03  9.469655e-03  9.208593e-02  2.977150e-01  6.957031e-01
# [16]  1.544891e+00  7.220152e+00  2.266144e+01  1.865402e+13  9.217375e+30
# [21]  5.552606e+45  1.105510e+47
# Warning: NAs produced in function "sqrt"
# >curvatureARMA
#       P   N Sign. level RMS Parameter RMS Intrinsic c*sqrt(F) Parameter
# [1,] 22 200        0.05  1.624603e+48  4.083931e+48        2.056528e+48
#      c*sqrt(F) Intrinsic Min Axis Ratio Max Axis Ratio
# [1,]        5.169704e+48       1.001183             NA

#  if(is.R())     tst <- c(22, 200, 0.05, 4.887000172697958e+24,
#      803263210.2793407,   6.1862808878576737e+24,  1016822523.0339906, 1, NaN)
#  if(is.Splus()) tst <- c(22, 200, 0.05, 2.19321390095902924, 
#      1.92482552083342284, 2.77631200307412485,     2.43656863335422758,1, NA)

#  if (print.values)
#     {cat("curvatureARMA\n") ; print.test.value(curvatureARMA, digits=digits)}
# These are very sensitive to simulation values ??
# R pre 1.0  tst <- c(22, 200, 0.05, 11316855942420566, 125981.7924284311,
#     14325608175410602,  159475.9007933074, 1)
#  if(is.S())  tst <- c(22, 200, 0.05,  12392476087969646,      155455.0076196983, 
#     15687197721934342,       196785.0027778072, 1)
  if(is.S())  tst <- c(22, 200, 0.05,  3.27852913374396824e+24,
    8.75713437411501169e+08, 4.15017421805572913e+24, 1.10853470629368877e+09,1)
  if(is.R()) tst <- c(22, 200, 0.05, 3.28476813494481106e+24,    26746456026841720,
      4.15807194928636139e+24,    33857393879614048, 1)

  all.ok <- test(curvatureARMA[-9], tst, all.ok, flag="15", fuzz=fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (synopsis) 
    {if (verbose) cat("All curvature tests completed")
     if (all.ok$ok) cat(" OK\n")
     else    
       {cat(", some FAILED!")
        if((!is.na(all.ok$error)) && (all.ok$error > fuzz.small))
            cat(" ( max. error magnitude= ", all.ok$error,")")
        cat("\n")
       }
    }
  invisible(all.ok)
}



dsecurvature.function.testsB <- function( verbose=T, synopsis=T, fuzz.small=1e-12, 
    fuzz.large=1e-6, show.details=F, show.extra.details=F,
    print.values=T, digits=18)
{# Tests of DSE curvature functions
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))


  all.ok <- list(ok=T, error=NA)

test <- function(tst, good, all.ok, flag="", fuzz=1e-14, print.values=F)
  {# initialize all.ok <- list(ok=T, error=NA)
   error <- max(abs(good-tst))
   if (fuzz < error) 
     {all.ok$error <- if(is.na(all.ok$error)) error else max(error, all.ok$error)
      all.ok$ok <- F  # assumes ok is initialized T for this to work 
      cat(flag, "failed! error = ", error,"\n")
      if (print.values) print.test.value(c(tst), digits=18)
     }
   all.ok
  }

  if (synopsis & !verbose) cat("All curvature tests ...")

#if(value < R_NSize || value > 5000000) changed from 1M in unix/system.c

# data size affects memory constraints
  data <- eg1.DSE.data.diff
   input.data(data) <- NULL
  output.data(data) <- output.data(data)[1:50,1:2]

  if (verbose) cat("DSE curvature test B 1 ...")
  VARmodel <- est.VARX.ls(data, re.add.means=F)
  SSmodel  <- l(to.SS(VARmodel),  data)
  ARMAmodel<- l(to.ARMA(SSmodel), data)

  spanVAR <- span(VARmodel)  
  tst <- c(0.114449869113756347, 0.114449869113496860, 0.073392835320007246,
      0.073392835319901109, 0.066734407677220553, 0.066734407677178142,
      0.066614701447440999, 0.066614701447334335, 0.063011563098793522,
      0.063011563098589879, 0.059329122705575062, 0.059329122705415301,
      0.056826595787831717, 0.056826595787806411, 0.046204330863178059,
      0.046204330863014273, 0.040188560079512797, 0.040188560079447523,
      0.034056955674954058, 0.034056955674626722, 0.028893589372264245,
      0.028893589372176717, 0.025983005223889789, 0.025983005223849501  )
  all.ok <- test(spanVAR, tst, all.ok, flag="1", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 2 ...")
  spanSS <- span(SSmodel)
  tst <-  
    c(0.10274860855471462662, 0.10274860855462775167, 0.06007019888802080099
     , 0.06007019888791977069, 0.06001503393765705852, 0.06001503393755355797
     , 0.05967088443568880490, 0.05967088443514689811, 0.05668932601611959693
     , 0.05668932601604840388, 0.05452098535006656699, 0.05452098534998019858
     , 0.05124584470731020913, 0.05124584470707382877, 0.05113559933370873806
     , 0.05113559933361908061, 0.04955894474981626524, 0.04955894474939047389
     , 0.04729297332666777126, 0.04729297332610608862, 0.04671869698239693863
     , 0.04671869698201638194, 0.04454781954178882453, 0.04454781954164353103
     , 0.04048455062849974639, 0.04048455062837719859, 0.02942866622217246361
     , 0.02942866622200048313, 0.01693320051034891832, 0.01693320051021537237
     , 0.01308641053414394059, 0.01308641053403496353, 0.01181519516958567939
     , 0.01181519516908161559, 0.00541747854002755207, 0.00541747853967942689
     , 0.00529351434412992337, 0.00529351434394933953, 0.00491893097942279044
     , 0.00491893097926101534, 0.00342671831624704161, 0.00342671831567546801
     , 0.00227476744198765396, 0.00227476744187563897, 0.00012162189247796396
     , 0.00012162189207734090, 0.00003274568362696781, 0.00003274568361473657)

  all.ok <- test(spanSS, tst, all.ok, flag="2", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 3 ...")
  spanARMA <- span(ARMAmodel)
  if (is.R()) tst <-
  c(4.9049410781519936e-01, 1.2425030501882078e-01, 1.1246172533887255e-01,
    1.0529465638718155e-01, 9.9992325730457826e-02, 9.3539845628889493e-02,
    9.2174599397343862e-02, 8.9496247916628913e-02, 8.4586129130414314e-02,
    8.0765436606922522e-02, 7.8623725634452205e-02, 7.4094697970647883e-02,
    6.8118147328717266e-02, 6.6571019509098564e-02, 6.5347673219782007e-02,
    6.2868363407609373e-02, 6.0029353951002326e-02, 5.8265043911548742e-02,
    5.7428063879464482e-02, 5.5270247854539004e-02, 5.3574863616695417e-02,
    4.9500699356060406e-02, 4.7721177486595059e-02, 4.5675211105535969e-02,
    4.5187209000554139e-02, 4.4123164386888028e-02, 4.2609744901397094e-02,
    4.0728262066476652e-02, 3.8500749322484873e-02, 3.6752466902023120e-02,
    3.6020724107593854e-02, 3.4123385578064033e-02, 3.3321675360527339e-02,
    2.9526045053003675e-02, 2.4664231456454182e-02, 2.4535892553740390e-02,
    1.9011587409302790e-02, 1.0185539990175529e-02, 7.9577476414227771e-03,
    4.4117567759992945e-03, 3.5579116575524851e-03, 2.5973873441597017e-03,
    2.3124010577862367e-03, 1.6333230757236616e-03, 1.3482756318466882e-03,
    4.0276257101139824e-04, 9.1217773859851665e-05, 8.4773331050698041e-05,
    7.5783716183056057e-09, 4.5904401626724121e-10, 1.1070461766649269e-17,
    9.4242715041701337e-18, 8.2559055659668990e-18, 7.6624154284980461e-18,
    6.8081425995733192e-18, 4.0398710160059630e-18, 6.6123246315825315e-19,
    3.8189436573397887e-19, 1.5940242413042422e-19, 9.4133611638589657e-20,
    7.4869852566030682e-20, 6.0095202032585055e-34, 4.7037657003150257e-34,
    3.5938593583744286e-34, 2.8591054742395741e-35, 0.0000000000000000e+00,
    0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
    0.0000000000000000e+00, 0.0000000000000000e+00 )

  if (is.Splus()) tst <-
    c(1.0034047190656252e+00, 4.8850479622289339e-01, 1.2408559806984479e-01,
     ,1.1241084023931164e-01, 1.0524964938676544e-01, 9.5268480653867810e-02,
     ,9.2759311159151175e-02, 9.1819598877975070e-02, 8.9050273189422796e-02,
     ,8.3629646994090756e-02, 8.0627890823442142e-02, 7.7853771461977417e-02,
     ,7.3417058329559359e-02, 6.8115229181312409e-02, 6.6487574775202618e-02,
     ,6.5268877320737190e-02, 6.2804879512010081e-02, 5.9798663205619551e-02,
     ,5.8138919962726084e-02, 5.7002371757449002e-02, 5.4146518855642631e-02,
     ,5.2602825342313142e-02, 4.9441695936262117e-02, 4.7269740278206025e-02,
     ,4.5517604797007144e-02, 4.5143602143697413e-02, 4.3842221290001454e-02,
     ,4.0984110441742688e-02, 3.9989100730095282e-02, 3.8350089480969470e-02,
     ,3.6071732859165054e-02, 3.5742661867439247e-02, 3.3945300390401167e-02,
     ,3.3080230063434836e-02, 2.8791495436511692e-02, 2.4595520370852313e-02,
     ,2.2520600570674183e-02, 1.8727379899113081e-02, 9.8908734880565584e-03,
     ,7.7946301146879322e-03, 4.3757311212158376e-03, 3.5514904518346609e-03,
     ,2.5730688907734350e-03, 2.2937687921734265e-03, 1.6334339201785643e-03,
     ,1.3490766143686318e-03, 5.2025361050217576e-04, 3.9749496568613595e-04,
     ,8.7754565756829609e-05, 8.4768386053768149e-05, 1.9920899452363507e-09,
     ,1.3537769691006131e-17, 1.2204712468141044e-17, 8.9797330692404288e-18,
     ,7.7118704229531706e-18, 6.6347042930705941e-18, 5.6925127578653899e-18,
     ,1.7886910482621617e-18, 1.1187688280849653e-18, 2.2384834985967148e-19,
     ,1.7676728214211909e-19, 1.4038899180867549e-19, 8.0812882075936788e-20,
     ,9.2192393808022497e-34, 6.6365234946111944e-34, 4.2858059029140618e-34,
     ,0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
     ,0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00 )

# Not clear why these should be different but now:
  tst <- c(  0.490494752810582801,  0.124254301561514971,  0.112464007925927312,
    0.105296151523361983,  0.0999930106260916302,  0.0935415658736187122,  
    0.0921756236732814843,  0.0894980196142749346,  0.084587879405886357,  
    0.080767820977991206,  0.0786269091436579692,  0.0740982050811652415,  
    0.0681193916577149394,  0.0665714610141954349,  0.0653501433331552389,  
    0.0628846530224142358,  0.0600354127079369554,  0.0582708846001192896,  
    0.0574305754952406169,  0.0552710746500495523,  0.053574894712188792,  
    0.0495021303437292587,  0.0477217132651450854,  0.0456814513508665004,  
    0.0451878212803662471,  0.0441241488416704852,  0.042609890004122658,  
    0.0407357969126947092,  0.0385011094551413538,  0.0367546517958144414,  
    0.0360356097631223324,  0.0341237444130534223,  0.033322353549055661,  
    0.029531589165205166,  0.024682000709630933,  0.0245390284959757336,  
    0.0190136770995758299,  0.010294991264683126,  0.00801491074106989478,  
    0.0045410856260552435,  0.00357647859566593597,  0.00308220555140791344,  
    0.00296860623932166312,  0.00248906796146492609,  0.00213323534758624389,  
    0.00157837322687163963,  0.00126529950976378707,  0.000400334715145694473,
    9.06116047208433366e-05,  8.46747193312290318e-05,  5.63955802833271451e-06,
    7.06607834857295999e-07,  1.02427642235884008e-17,  8.2945253576028173e-18,
    7.73322892221025617e-18,  5.87799494460792363e-18,  5.34659049636090042e-18,
    2.16262032052361535e-18,  1.36681894882810078e-18,  2.42341824448017461e-19,
    2.37700839985245846e-19,  1.85660348662921073e-19,  1.83180499055525642e-19,
    6.63065315055060987e-34,  6.02668158348751715e-34,  3.59344821900416301e-34,
    0,                    0,                    0,                    0, 0)
  all.ok <- test(spanARMA, tst, all.ok, flag="3", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 4 ...")
  curvatureVAR <- curvature(VARmodel)$stats
  all.ok <- test(curvatureVAR, c(24, 100, 0.05, 0, 0, 0, 0, 1, 1 ), all.ok, flag="4", fuzz=1e-5,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 5 ...")
  curvatureSS <- curvature(SSmodel, warn=F)$stats
#  if (print.values) print.test.value(curvatureSS, digits=digits)

# if(is.R())     tst <- c(48, 100, 0.05, 323.99227471499745, 124.74717834454975,
#        409.3130651286905,  157.59835625487983, 1.0000000021695887, NaN)[-9]
# if(is.Splus()) tst <- c(48, 100, 0.05, 255.387386434704,    98.33696083718515,
#        322.6416248003061, 124.23321787877502, 1, 1.0000000004520695  ) [-9]

  tst <-c(48, 100, 0.0500000000000000028, 323.992363637504695,
          124.759997096571936,  409.313177468233278,  157.614550723361958,
            1.0000000014689332)
#  ok <- fuzz.small >  max(abs( curvatureSS[-9] -  tst))

#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(curvatureSS[-9], tst, all.ok, flag="5", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 6 ...")
  curvatureARMA <- curvature(ARMAmodel, warn=F)$stats
#  if (print.values) print.test.value(curvatureARMA, digits=digits)
# if(is.R())     tst <- c(71,100,0.05, 8.083857768891898e+23, 354675866.7879350,
#          1.0653907038344889e+24, 467435699.82689363, 1, NaN )[-9]
# if(is.Splus()) tst <- c(72,100,0.05, 31947799733885313024, 1708051.5249938569,
#       42274456907383595008,  2260154.101077503, 1, 1.0000267463352852 )[-9] 

# tests looks suspicious

  tst <- c(71, 100, 0.0500000000000000028, 1.48023299583791679e+23, 
           7.16541930125915901e+21,  1.95083401806432887e+23,  
           9.443475294697275e+21,   1 )

#  ok <- fuzz.small >  max(abs( curvatureARMA[-9] -  tst))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(curvatureARMA[-9]*c(1,1,1,1e-23, 1e-21, 1e-23, 1e-21, 1),
                           tst*c(1,1,1,1e-23, 1e-21, 1e-23, 1e-21, 1),
		 all.ok, flag="6", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 7 ...")
  hessianVAR <- hessian(VARmodel)
#  ok <- fuzz.small >  max(abs( sum(hessianVAR) - 7219.22210394543526 ))
#            if(is.Splus()) 7219.717083137912 else if(is.R()) 7219.19366223377))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(sum(hessianVAR), 7219.22210394543526, all.ok, flag="7", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 8 ...")
  hessianSS <- hessian(SSmodel)

# previously if(is.Splus()) 7844.3395239153897 else if(is.R()) 7841.271986002))

  all.ok <- test(sum(hessianSS), 7841.24813340843411, all.ok, flag="8", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test B 9 ...")
  hessianARMA <- hessian(ARMAmodel)
#  if (print.values) print.test.value(sum(hessianARMA), digits=digits)

#  ok <- fuzz.small >  max(abs( sum(hessianARMA) -  256440.198630697385  ))
  #  +  if(is.Splus()) 1789846677.6677122 else if(is.R()) 90636.84015934517))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(sum(hessianARMA), 256440.198630697385, all.ok, flag="9", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (synopsis) 
    {if (verbose) cat("All curvature tests completed")
     if (all.ok$ok) cat(" OK\n")
     else    
       {cat(", some FAILED!")
        if((!is.na(all.ok$error)) && (all.ok$error > fuzz.small))
            cat(" ( max. error magnitude= ", all.ok$error,")")
        cat("\n")
       }
    }
  invisible(all.ok)
}


# following from dse1 tests are much slower

dsecurvature.function.testsC <- function( verbose=T, synopsis=T, fuzz.small=1e-12, 
    fuzz.large=1e-6, show.details=F, show.extra.details=F,
    print.values=T, digits=18)
{# Tests of DSE curvature functions
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  all.ok <- list(ok=T, error=NA)

test <- function(tst, good, all.ok, flag="", fuzz=1e-14, print.values=F)
  {# initialize all.ok <- list(ok=T, error=NA)
   error <- max(abs(good-tst))
   if (fuzz < error) 
     {all.ok$error <- if(is.na(all.ok$error)) error else max(error, all.ok$error)
      all.ok$ok <- F  # assumes ok is initialized T for this to work 
      cat(flag, "failed! error = ", error,"\n")
      if (print.values) print.test.value(c(tst), digits=18)
     }
   all.ok
  }


  if (synopsis & !verbose) cat("All curvature tests ...")

  if (verbose) cat("DSE curvature tests C...\n")
  if (verbose) cat("DSE curvature test C 1 ...")
  VARmodel <- est.VARX.ls(eg1.DSE.data.diff, re.add.means=F)
  SSmodel  <- l(to.SS(VARmodel),  eg1.DSE.data.diff)
  ARMAmodel<- l(to.ARMA(SSmodel), eg1.DSE.data.diff)

  spanVAR <- span(VARmodel)

#  if (print.values) print.test.value(spanVAR, digits=digits)
  # either S or R values work here with fuzz.large but try with something tighter
  # R 1.0.1:
#  if (is.R())
  tst <- c(
     431.324809764142003,  431.324809764133477,  431.324809764125973,  
     29.9938661246376412,  29.9938661245786164,  29.9938661245442368,  
     14.0487108329584771,   14.048710832937374,  14.0487108329158925, 
      8.41716993883251874,  8.41716993883184195,  8.41716993881048836,   
      6.6727488248550193,  6.67274882470962183,  6.67274882470017783,  
      6.23686447689305545,  6.23686447689176671,  6.23686447680381217,  
      0.278224057276565351,  0.278224057275740233,  0.278224057275401837,  
      0.226296510339683066,  0.226296510339373452,  0.226296510338978879,  
      0.22109512238799589,  0.221095122387542337,  0.221095122387375581,  
      0.217635227413440741,  0.217635227413369131,  0.217635227412733029,  
      0.210139509606726088,  0.210139509606644848,  0.210139509606350583,  
      0.193542574580348015,  0.193542574580293447,  0.193542574579900289,  
      0.133100015736768662,  0.133100015736627497,  0.133100015736183325,  
      0.124513362052187021,  0.124513362052168133,  0.124513362051330567,  
      0.121954031479664826,  0.121954031479572386,  0.121954031479351521,  
      0.113732578253665478,  0.113732578253511143,  0.113732578253048111,  
      0.089473111199309549,  0.0894731111990460237,  0.0894731111990405698,  
      0.0871822427356184065,  0.0871822427354376067,  0.0871822427353293739,  
      0.0777910835236551579,  0.0777910835236260007,  0.0777910835233898562,  
      0.0595266223042596881,  0.0595266223042083611,  0.0595266223040415363,  
      0.0499760073414238712,  0.0499760073413548916,  0.0499760073413023087,  
      0.0485695297799133355,  0.0485695297798213396,  0.048569529779736699,  
      0.0474465370323315233,  0.0474465370321883045,  0.0474465370320783369,  
      0.0443198146313130703,  0.0443198146313004693,  0.0443198146312605984 )

#  ok <- fuzz.small >  max(abs( spanVAR - tst ))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(spanVAR, tst, all.ok, flag="1", fuzz=100*fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 2 ...")

  spanSS <- span(SSmodel)
#  if (print.values) print.test.value(spanSS, digits=digits)
#  ok <- ok & fuzz.small >  max(abs( spanSS -  tst))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(spanSS, spanVAR, all.ok, flag="2", fuzz=10*fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 3 ...")
  spanARMA <- span(ARMAmodel)
  all.ok <- test(spanARMA, 0, all.ok, flag="3", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 4 ...")
  curvatureVAR <- curvature(VARmodel)$stats
  tst <- c(24, 100, 0.05, 1.4877780902309642e-06, 3.2049317816656057e-06, 
           1.917679915995849e-06, 4.1310147999845919e-06, 1, 1           )
  all.ok <- test(curvatureVAR, tst, all.ok, flag="4", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 5 ...")
  curvatureSS <- curvature(SSmodel)$stats
  tst <- c(48, 100, 0.05, 255.387386434704, 98.33696083718515,
           322.6416248003061, 124.23321787877502, 1, 1.0000000004520695  )
  all.ok <- test(curvatureSS, tst, all.ok, flag="5", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 6 ...")
  curvatureARMA <- curvature(ARMAmodel)$stats
  tst <-c(72, 100, 0.05, 31947799733885313024, 1708051.5249938569,
   42274456907383595008,  2260154.101077503, 1, 1.0000267463352852 )
# above looks suspicious
  all.ok <- test(curvatureARMA, tst, all.ok, flag="6", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 7 ...")
  hessianVAR <- hessian(VARmodel)
  all.ok <- test(sum(hessianVAR), 7219.717083137912, all.ok, flag="7", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 8 ...")
  hessianSS <- hessian(SSmodel)
  all.ok <- test(sum(hessianSS), 7844.3395239153897, all.ok, flag="8", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 9 ...")
  hessianARMA <- hessian(ARMAmodel)
  all.ok <- test(sum(hessianARMA), 1789846677.6677122, all.ok, flag="9", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (synopsis) 
    {if (verbose) cat("All curvature tests completed")
     if (all.ok$ok) cat(" OK\n")
     else cat(", some FAILED! all.ok$error = ", all.ok$error,"\n")
    }
  invisible(all.ok)
}

#######################################################################

#                    end

#######################################################################


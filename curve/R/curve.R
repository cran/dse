#   2000/04/18 12:08:12 


############################################################################

#    functions for gradient calculation

############################################################################

numerical.grad <- function(func,x, eps=1e-12) {
#  very simple (crude) numerical approximation
  f <-func(x)
  df <-1:length(x)
  for (i in 1:length(x)) {
    dx <- x
    dx[i] <- dx[i] +eps 
    df[i] <- (func(dx)-f)/eps
   }
df
}


richardson.grad <- function(func, x, d=0.01, eps=1e-4, r=6, show.details=F){
# This function calculates a numerical approximation of the first
#   derivative of func at the point x. The calculation
#   is done by Richardson's extrapolation (see eg. G.R.Linfield and J.E.T.Penny
#   "Microcomputers in Numerical Analysis"). The method should be used if
#   accuracy, as opposed to speed, is important.
#
#  *  modified by Paul Gilbert from orginal code by XINGQIAO LIU.
# CALCULATES THE FIRST ORDER DERIVATIVE 
#     VECTOR using a Richardson  extrapolation.
#
#  GENERAL APPROACH
#     --  GIVEN THE FOLLOWING INITIAL VALUES:
#             INTERVAL VALUE D, NUMBER OF ITERATIONS R, AND
#             REDUCED FACTOR V.
#      - THE FIRST ORDER aproximation to the DERIVATIVE WITH RESPECT TO Xi IS
#
#           F'(Xi)={F(X1,...,Xi+D,...,Xn) - F(X1,...,Xi-D,...,Xn)}/(2*D)
#       
#     --  REPEAT r TIMES  with successively smaller D  and 
#          then apply Richardson extraplolation
#
#  INPUT
#       func    Name of the function.
#       x       The parameters of func.
#       d       Initial interval value (real) by default set to 0.01*x or
#               eps if x is 0.0.
#       r       The number of Richardson improvement iterations.
#       show.details    If T show intermediate results.
#  OUTPUT
#
#       The gradient vector.

  v <- 2               # reduction factor.
  n <- length(x)       # Integer, number of variables.
  a.mtr <- matrix(1, r, n) 
  b.mtr <- matrix(1, (r - 1), n)
#------------------------------------------------------------------------
# 1 Calculate the derivative formula given in 'GENERAL APPROACH' section.
#   --  The first order derivatives are stored in the matrix a.mtr[k,i], 
#        where the indexing variables k for rows(1 to r),  i for columns
#       (1 to n),  r is the number of iterations, and n is the number of
#       variables.
#-------------------------------------------------------------------------  

  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  { # successively reduce h                
     for(i in 1:n)  {
         x1.vct <- x2.vct <- x
         x1.vct[i]  <- x[i] + h[i]
         x2.vct[i]  <- x[i] - h[i]
         if(k == 1) a.mtr[k,i] <- (func(x1.vct) - func(x2.vct))/(2*h[i])
         else{
           if(abs(a.mtr[(k-1),i])>1e-20)
                 # some functions are unstable near 0.0              
                 a.mtr[k,i] <- (func(x1.vct)-func(x2.vct))/(2*h[i])
            else  a.mtr[k, i] <- 0
          }
      }
     h <- h/v     # Reduced h by 1/v.
    }	
   if(show.details)  {
        cat("\n","first order approximations", "\n")		
        print(a.mtr, 12)
    }

#------------------------------------------------------------------------
# 1 Applying Richardson Extrapolation to improve the accuracy of 
#   the first and second order derivatives. The algorithm as follows:
#
#   --  For each column of the 1st and 2nd order derivatives matrix a.mtr,
#       say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
#       new sequence of approximations B1, B2, ..., Br used the formula
#
#          B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
#
#             N.B. This formula assumes v=2.
#
#   -- Initially m is taken as 1  and then the process is repeated 
#      restarting with the latest improved values and increasing the 
#      value of m by one each until m equals r-1
#
# 2 Display the improved derivatives for each
#   m from 1 to r-1 if the argument show.details=T.
#
# 3 Return the final improved  derivative vector.
#-------------------------------------------------------------------------

  for(m in 1:(r - 1)) {		
#     for(i in 1:(r - m)) b.mtr[i,]<- (a.mtr[(i+1),]*(4^m)-a.mtr[i,])/(4^m-1)
#     a.mtr<- b.mtr
     a.mtr<- (a.mtr[2:(r+1-m),,drop=F]*(4^m)-a.mtr[1:(r-m),,drop=F])/(4^m-1)
     if(show.details & m!=(r-1) )  {
        cat("\n","Richarson improvement group No. ", m, "\n")		
        print(a.mtr[1:(r-m),,drop=F], 12)
      }
   }
a.mtr
}


############################################################################

#    functions for projecting one space on another

############################################################################



project <- function(c1, c2, signif = 0.05, eps=1e-5, fuzz=1e-14,
    show.details=F, warn=T)
{  p1 <- c1$Dlist$p
   N <- length(c1$Dlist$f0)   # dimesion of sample space
   p2 <- c2$Dlist$p
   if (N != length(c2$Dlist$f0)) stop("sample space dimesions do not agree.")
   if (p2 < p1) stop("submodel should be the first argument.")
   if (max(abs(c1$Dlist$f0 - c2$Dlist$f0)) > fuzz)
      warning("projection not computed at the same sample space point as the original curvature calculation.")

   residual <- c(c1$Dlist$f0)   # c() in case it is a column vector
   s.sqr <- sum(residual^2)/(N-2)  # check for time series??

   cat("p1=",p1," p2=",p2,"N=",N,"\n")
   cat("max. span of acceleration space as determined by number of vectors \nin 2nd derivative array:    ")
   cat("   1st model=", p1 *(1+p1)/2,
       ",2nd model=", p2 *(1+p2)/2,"\n" )   

   d1 <- svd(c1$Dlist$D[,1:p1])$d
   d2 <- svd(c2$Dlist$D[,1:p2])$d
   dj <- svd(cbind(c1$Dlist$D[,1:p1], c2$Dlist$D[,1:p2]))$d
   cat("span of tangent vectors:       1st model=", sum(d1 > eps*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*d2[1]))
   cat(", jointly=",       sum(dj > eps*dj[1]), "\n")
 
   d1 <- svd(c1$Dlist$D[,-(1:p1)])$d
   d2 <- svd(c2$Dlist$D[,-(1:p2)])$d
   dj <- svd(cbind(c1$Dlist$D[,-(1:p1)], c2$Dlist$D[,-(1:p2)]))$d
   cat("span of acceleration vectors:  1st model=", sum(d1 > eps*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*d2[1]))
   cat(", jointly=",       sum(dj > eps*dj[1]), "\n")
   cat("using eps*100 (sensitivity for significance of svd$d )\n")
   cat("span of acceleration vectors:  1st model=", sum(d1 > eps*100*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*100*d2[1]))
   cat(", jointly=",       sum(dj > eps*100*dj[1]), "\n")

   m2 <- p2 * (p2 + 3)/2  # m= p+pp is second dimension of D.
   m2 <- min(m2,dim(c2$Dlist$D)[2]) 
   m1 <- p1 * (p1 + 3)/2 

#  c1 projected into c2

   QRofD2 <- qr(c2$Dlist$D)
   T1inT2 <- qr.qty(QRofD2, c1$Dlist$D[,  1:p1 ] )[  1:p2 ,,drop=F]  
   N1inT2 <- qr.qty(QRofD2, c1$Dlist$D[,-(1:p1)] )[  1:p2 ,,drop=F]  
   N1inN2 <- qr.qty(QRofD2, c1$Dlist$D[,-(1:p1)] )[(p2+1):m2,,drop=F]  

#  c2 projected onto c1

   QRofD1 <- qr(c1$Dlist$D)
   N2inN1      <- qr.qty(QRofD1, c2$Dlist$D[,-(1:p2)] )[(p1+1):m1,,drop=F]  
   T2inN1      <- qr.qty(QRofD1, c2$Dlist$D[,  1:p2 ] )[(p1+1):m1,,drop=F] 
   N2andT2inN1 <- qr.qty(QRofD1, c2$Dlist$D           )[(p1+1):m1,,drop=F]
   N2andT2inc1 <- qr.qty(QRofD1, c2$Dlist$D           )

   N1inN1 <- qr.qty(QRofD1, c1$Dlist$D[,(p1+1):m1] )[(p1+1):m1,,drop=F]  
#   v1 <- svd(N2andT2inT1)
#   v2 <- svd(T1inT1)
browser()

   cur2in1 <- rel.curvature(s.sqr,N2andT2inc1[1:p2,1:p1], N2andT2inc1[,(p1+1):m1], 
                      show.extra.details=show.details)
   C2in1 <- list(C.parameter=cur2in1[1:p1,,     ,drop=F],
                 C.intrinsic=cur2in1[(p1+1):m1,,,drop=F])

#    z <- rel.curvature(s.sqr,N1inN1[1:p1,1:p1], N2andT2inc1[,(p1+1):m1], 
#                       show.extra.details=show.details)
#    zz <- list(C.parameter=z[1:p1,,     ,drop=F],
#                  C.intrinsic=z[(p1+1):m1,,,drop=F])

# svd(   zz$C.intrinsic[1,,])$d
# svd( c1$C$C.intrinsic[1,,])$d
# svd( c2$C$C.intrinsic[1,,])$d
# svd(C2in1$C.intrinsic[1,,])$d

# svd( c1$C$C.intrinsic[1,,])$d / svd(C2in1$C.intrinsic[1,,])$d[1:2]

# eigen( c1$C$C.intrinsic[1,,])$values
# eigen( c2$C$C.intrinsic[1,,])$values
# eigen(C2in1$C.intrinsic[1,,])$values

   effective<-effective.curvature(cur1on2,QRofD2, residual, s.sqr,
                      show.details=show.details, warn=warn)

   cstats1on2 <-curvature.stats(cur1on2, N, signif=signif)

   R1 <- qr.qty(QRofD1, c2$Dlist$D )[1:m1,,drop=F]  
   cur1 <- rel.curvature(s.sqr,R1[1:p1,1:p1], R1[,(p1 + 1):m1], 
                      show.extra.details=show.details)
   effective1<-effective.curvature(cur1,QRofD1, residual, s.sqr,
                      show.details=show.details)

   cstats1 <-curvature.stats(cur1, N, signif=signif)
browser()

   result <- rbind(c1$stats,c(p2, N, signif, cstats1on2,
                   effective$min.axis.ratio, effective$max.axis.ratio),
                   c2$stats)
   dimnames(result)<-list(c("original c1", "projected", "original c2"),
                          c("Parms","Sample","Sign. level","RMS Parameter", 
         "RMS Intrinsic","c*sqrt(F) Parameter","c*sqrt(F) Intrinsic",
         "Min Axis Ratio","Max Axis Ratio")) 
   result
}


#######################################################################

#            calculate Hessian 

#######################################################################



hessian <- function (obj, ...)  UseMethod("hessian")


hessian.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4,r=6)
{  D <- genD.default(func,x,func.args=func.args, d=d, eps=eps, r=r)$D
   H <- diag(0,length(x))
   u <- 0
   for(i in 1:length(x))
     {for(j in 1:i) 
        {u <- u + 1
         H[i,j] <- D[,u]
     }  }
   H <- H +t(H)
   diag(H) <- diag(H)/2
   H
}


#######################################################################

#              span (dimension of tangent space) calculation

#######################################################################

span.v00 <- function(func, x,d=0.01, eps=1e-4,r=6, show.details=F){
#  Calculate tangent vectors at the point x. 
#   Note that func can be vector valued.
#
#   This function uses Richardson extrapolation  (for more details      
#      see the functions richardson.grad and genD) to get a numerical 
#      approximation of the tangent vectors to the parameter 
#      manifold. SVD is then used to calculate their span.     
#
#       func    the name of a function which returns the
#                residual vector for a given parameter vector.
#       x       The parameters of func.
#
#       d       initial interval value by
#               default set to 0.01*x or eps if x is 0.0 .
#       r       the number of repetions with successly smaller d. 
#       show.details    if TRUE then detailed calculations are shown.
#       v       reduction factor. This could be a parameter but ...
#               N.B. must be 2 for richarson formula as used.

  v <- 2
  p <- length(x)     #  number of parameters (= dim of M if not over parameterized.
  n <- len(func(x))  #  dimension of sample space.
  D <- array(1, c(n, p, r)) 
  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  # successively reduce h 
    {for(i in 1:p)  # each parameter  - first deriv.
       D[,i,k]<-(func(x+(i==(1:p))*h)-func(x-(i==(1:p))*h))/(2*h[i]) # F'(i)
     h <- h/v     # Reduced h by 1/v.
    }	
  if(show.details)  
      {cat("\n","first order approximations", "\n")		
       for(i in 1:p){ cat(" parameter " ,i,"\n");print(D[,i,(1:r)], 12)}
      }
  for(m in 1:(r - 1)) 		
     {D<- (D[,,2:(r+1-m)]*(4^m)-D[,,1:(r-m)])/(4^m-1)
      if(show.details & m!=(r-1) )  
        {cat("\n","Richarson improvement group No. ", m, "\n")		
         for(i in 1:p){ cat(" parameter " ,i,"\n"); print(D[,i,(1:(r-m))], 12) }
        }
      }
   svd(D)$d
}



span.v0 <- function(func, x,d=0.01, eps=1e-4,r=6, show.details=F){
#  span performs a svd of the tangent vectors at the point x. 
#   (Note that func can be vector valued. This can be used
#   to calculate the dimension of the tangent space (ie. by over specifying 
#   the model and counting the number of significant singular values).
#  The singular values are returned.
#
#   This function uses Richardson extrapolation  (for more details      
#      see the functions richardson.grad and genD) to get a numerical 
#      approximation of the tangent vectors to the parameter 
#      manifold. SVD is then used to calculate their span.     
#
#       func    the name of a function which returns the
#                residual vector for a given parameter vector.
#       x       The parameters of func.
#
#       d       initial interval value by
#               default set to 0.01*x or eps if x is 0.0 .
#       r       the number of repetions with successly smaller d. 
#       show.details    if TRUE then detailed calculations are shown.
#       v       reduction factor. This could be a parameter but ...
#               N.B. must be 2 for richarson formula as used.

  v <- 2
  p <- length(x)     #  number of parameters (= dim of M if not over parameterized.
  n <- len(func(x))  #  dimension of sample space.
  D <- array(1, c(n, p, r)) 
  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  # successively reduce h 
    {for(i in 1:p)  # each parameter  - first deriv.
       D[,i,k]<-(func(x+(i==(1:p))*h)-func(x-(i==(1:p))*h))/(2*h[i]) # F'(i)
     h <- h/v     # Reduced h by 1/v.
    }	
  if(show.details)  
      {cat("\n","first order approximations", "\n")		
       for(i in 1:p){ cat(" parameter " ,i,"\n");print(D[,i,(1:r)], 12)}
      }
  for(m in 1:(r - 1)) 		
     {D<- (D[,,2:(r+1-m)]*(4^m)-D[,,1:(r-m)])/(4^m-1)
      if(show.details & m!=(r-1) )  
        {cat("\n","Richarson improvement group No. ", m, "\n")		
         for(i in 1:p){ cat(" parameter " ,i,"\n"); print(D[,i,(1:(r-m))], 12) }
        }
      }
   svd(D)$d
}

span <- function (obj, ...)  UseMethod("span")

span.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4,r=6, show.details=F)
{
#   performs a svd of the tangent vectors at the point x. This can be used
#   to calculate the dimension of the tangent space (ie. by over specifying 
#   the model and counting the number of significant singular values).
#  The singular values are returned.
#
#   This function uses Richardson extrapolation  (for more details      
#      see the functions richardson.grad and genD) to get a numerical 
#      approximation of the tangent vectors to the parameter 
#      manifold. SVD is then used to calculate their span.     
#
#  func     a charater string of the name of a function which returns the
#                residual vector for a given parameter vector.
#   x       The parameters of func with respect to which derivative is calculated.
#  func.args A list of any other arguments to the function
#
#   d       initial interval value by
#               default set to 0.01*x or eps if x is 0.0 .
#   r       the number of repetions with successly smaller d. 
#   show.details    if TRUE then detailed calculations are shown.
#   v       reduction factor. This could be a parameter but ...
#             N.B. must be 2 for richarson formula as used.

  v <- 2
  p <- length(x)     #  number of parameters (= dim of M if not over parameterized.
  n <- length(do.call(func,append(list(x), func.args)))  #dimension of sample space.
  D <- array(1, c(n, p, r)) 
  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  # successively reduce h 
    {if(show.details) cat("\n k=",k, " p: ")
     for(i in 1:p)  # each parameter  - first deriv.
      {if(show.details) cat(" ",i)
       D[,i,k]<-(do.call(func,append(list(x+(i==(1:p))*h), func.args)) -
                 do.call(func,append(list(x-(i==(1:p))*h), func.args)))/(2*h[i])
       NULL 
      }
     h <- h/v     # Reduced h by 1/v.
    NULL 
    }	
  if(show.details)  
      {cat("\n","first order approximations", "\n")		
       for(i in 1:p){ cat(" parameter " ,i,"\n");print(D[,i,(1:r)], 12); NULL}
      }
  for(m in 1:(r - 1)) 		
     {D<- (D[,,2:(r+1-m)]*(4^m)-D[,,1:(r-m)])/(4^m-1)
      if(show.details & m!=(r-1) )  
        {cat("\n","Richarson improvement group No. ", m, "\n")		
         for(i in 1:p){ cat(" parameter " ,i,"\n"); print(D[,i,(1:(r-m))], 12) }
        }
      NULL
     }
   svd(D)$d
}


#######################################################################

#                            curvature calculation

#######################################################################




print.curvature <- function(x, ...)  { print(x$stats, ...) }

curvature <- function (obj, ...)  
{ # calculate the Bates and Watts intrinsic and parameter effects curvature.
  UseMethod("curvature")
}


curvature.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4,r=6,
      signif=0.05, show.details=F, warn=T)
{curvature(genD.default(func,x, func.args=func.args, d=d, eps=eps,r=r),
     signif=0.05, show.details=F, warn=warn)
}


curvature.Darray <- function(Dlist, signif = 0.05, show.extra.details=F, show.details=show.extra.details, warn=T)
{
# Curvature summary statistics as in Bates and Watts.
#  Dlist is a list as generated by genD.default and genD.TSestModel,
#       with the 3 elements as follows:
#   $D is a matrix of first(gradients) and second order partial
#      derivatives organized in the same manner as Bates and 
#      Watts. (first p columns are the gradients and the 
#      next p(p+1)/2 columns are the lower triangle of the Hessian).
#   $p is the dimension of the parameter space=dim of the tangent space.
#   $f0 is the function value at the point where the matrix D 
#        was calculated. (The curvature calculation should not/does not? depend 
#        on this value - but it should be the right dimension and 0's do
#        not work.

#  Page references are to Bates and Watts(1983), Nonlinear Regression 
#   Analysis and Its Applications.
#   Examples are also from this text. e.g.- model A data set 3.

#   modified substantially by Paul Gilbert (May, 1992 and Dec, 1996)
#    from original code by Xingqiao Liu,   May, 1991.
#
   if (all(Dlist$D==0))
     stop("The array of first and second order derivatives is all zeros.")

   p <- Dlist$p
   n <- dim(Dlist$D)[1]   # dimesion of sample space
   pp <- p*(p+1)/2 # p' max number of independent col in the acceleration array.
   m <- p * (p + 3)/2  # m= p+pp is second dimension of D. It is also the max 
                       # span of the combined tangent and accelation spaces.
   if (m!=dim(Dlist$D)[2]) stop("Dimesion of D matrix is not consistent with p")
   residual <- c(Dlist$f0)   # c() in case it is a column vector
   s.sqr <- sum(residual^2)/(length(residual)-2)
   if(show.extra.details) cat("p=",p," pp=",pp," m=",m, " n=",n,"\n")
 
#  Form QR decomposition of D and produce matrix Q and R.
#

   QRofD <- qr(Dlist$D) 

# Now calculate R, which is D rotated onto an othogonal basis with the 
# tangent space first, then acceleration space. 
# In this basis first rows of R are in the tangent and the next rows are in the
# acceleration space. Following rows would be
# all zeros  and so those rows are truncated from R.

#   Q <- qr.qy(QRofD, diag(1,n)) [,1:m]
#   R <- t(Q) %*% Dlist$D   this is the same as
#  R[1:m,]could be used to truncate but it is possible for m to exceed dim(D)[1]
   R <- qr.qty(QRofD, Dlist$D )  

   if(show.extra.details) 
        {cat( "Matrix D, tangent and acceleration vectors in columns ")
         cat("organized as in Bates and Watts p235, 236\n")
         print(Dlist$D)
         cat("Residual Vector\n"); print(residual)
	 cat("Matrix Q from QR decomposition of D\n")
         print( qr.qy(QRofD, diag(n)))
	}
   if(show.details) 
        {cat("Matrix R (organized as in Bates and Watts R1 p236)\n")
         print(R)
	}

#  The matrix R partitions into p by p matrix R11 and m by p'  matrix R2.
#  R2 is Bates and Watts  R12 
#                         R22  p236

#   R11 <- R[1:p,1:p]
#   R2 <- R[1:m,(p + 1):m] 

#                                     R[1:m,(p + 1):m]
   cur <- rel.curvature(s.sqr,R[1:p,1:p], R[,(p + 1):m], 
                      show.extra.details=show.extra.details)

   if (m>dim(cur)[1]) 
     {m <- dim(cur)[1]
      warning("acceration space dimension reduced for smaller sample space.")
     }
   C <- classed(list(C.parameter=cur[1:p,,,drop=F],# curvatureArray constructor
             C.intrinsic=cur[(p+1):m,,,drop=F]),"curvatureArray" )

# effective residual curvature. Bates and Watts p260
# Calculate the scaled RMS curvatures relative to a confidence disk 
# radius and extreme axis ratios. ref. Bates and Watts p254 eqn (7.23).
#  and Bates and Watts J.R. Statist.Soc. B (1980).

   effective<-effective.curvature(cur,QRofD, residual, s.sqr, 
                              show.details=show.extra.details, warn=warn)
   cstats <-curvature.stats(cur, n, signif=signif)

   result <- matrix(c(p, n, signif, cstats,
               effective$min.axis.ratio, effective$max.axis.ratio),1,9)
   dimnames(result)<-list(NULL,c("Parms","Sample","Sign. level","RMS Parameter",
         "RMS Intrinsic","c*sqrt(F) Parameter","c*sqrt(F) Intrinsic",
         "Min Axis Ratio","Max Axis Ratio")) 
   classed(list(stats=result, Dlist=Dlist, C=C), "curvature")# curvature.Darray constructor
}

print.curvatureArray <- function(x, ...)
  {cat("Parameter effects relative curvature array:\n")
   for (i in 1:dim(x$C.parameter)[1]) 
       {print(x$C.parameter[i,,], ...); cat("\n"); NULL}
   cat("Intrinsic relative curvature array:\n")
   for (i in 1:dim(x$C.intrinsic)[1])
       {print(x$C.intrinsic[i,,], ...); cat("\n"); NULL}
   invisible(x)
  }



curvature.stats <- function(cur, n, signif=0.05)
{# See documentation in in curvature.Darray
 p  <- dim(cur)[2]
 pp <- dim(cur)[1] - p
 cstats <- rep(NA,4)
   #  i from 1 to p for the RMS parameter effects curvature.
   #  i from p+1 to m=p+p' for the RMS intrinsic effects curvature.
   di <- 0
   for(i in 1:p) di <- di + sum(diag(cur[i,,]))^2
   cstats[1] <-  2*sum(cur[1:p,,]^2) + di
   cstats[1] <-  sqrt(cstats[1]/(p*(p+2)))
   di <- 0
   for(i in p+(1:pp)) di <- di + sum(diag(cur[i,,]))^2
   cstats[2] <-  2*sum(cur[p+(1:pp),,]^2) + di
   cstats[2] <-  sqrt(cstats[2]/(p*(p+2)))
  value.f <- qf(1 - signif, p, n - p)	# F(p, n-p; signif)
  cstats[3] <- cstats[1] * sqrt(value.f) # c*sqrt(F) Parameter
  cstats[4] <- cstats[2] * sqrt(value.f) # c*sqrt(F) Intrinsic
  cstats
}


effective.curvature <- function(cur,QRofD, residual, s.sqr,
			 show.details=F, warn=T)
{# Transform the residual vector by multiply Q-transpose and sqrt(s.sqr*p).
 # Calculate the p by p effective residual curvature matrix B and its 
 #  eigenvalues. Bates and Watts p260
 p  <- dim(cur)[2]
 pp <- dim(cur)[1] - p
 Qz <- qr.qty(QRofD, residual)[p+(1:pp)]/c(sqrt(p * s.sqr))  # p' vector of rotated residuals
 #pxp effective residual (intrinsic) curvature matrix. Bates and Watts p260
   # Qz is a vector of length pp.  
   B <- matrix(NA,p,p) 
   A.intrin <- cur[p+(1:pp),,]
 #   for (i in 1:p)
 #       for (j in 1:p)
 #         B[i,j] <-  A.intrin[,i,j] %*% Qz 
 # following is equivalent to above
   B <- apply(A.intrin * Qz, c(2,3), sum)
   if(show.details) 
     { cat("Scaled and Q-transformed Residual Vector\n"); print(Qz)
       cat("Effective residual curvature matrix B (Bates and Watts p 260)\n")
       print(B)
     }
   if ( max(abs(B- t(B))) > 1e-7 ) warning("B is not symmetric.")  #Rbug
   eigv <- eigen(B)$values
   if(max(Im(eigv)) >= 1e-15)
       warning("Calculated eigenvalues have imaginary part.")

   if(show.details)   {cat("Eigenvalues of B\n"); print(eigv )}

   if(warn && (max(Re(eigv)) >= 1.0))
     {warning( paste(
      "N.B. (I-B) is not positive definite as it should be at a local min! ",
      "Axis ratio will evaluate to NA or NaN. ",
      "Eigenvalues of B ", eigv))
     }

 #Mod in the following gets rid of small Im part, but the eigv should be real
# axis.ratios <- 1/sqrt(1-eigv)
# list(B=B, eigenvalues=eigv,
#      min.axis.ratio=min(axis.ratios),
#      max.axis.ratio=max(axis.ratios)  )
# Bates and Watts seem to use the following but the above should be considered.
# It seems possible to have an eigv < -1
# sqrt gives NaN for neg. numbers but test avoids producing additional warnings.
 list(B=B, eigenvalues=eigv,
      min.axis.ratio= if(1 > min(Mod(eigv))) 1/sqrt(1-min(Mod(eigv))) else NaN,
      max.axis.ratio= if(1 > max(Mod(eigv))) 1/sqrt(1-max(Mod(eigv))) else NaN)
}

rel.curvature <- function(s.sqr,R11,R2, show.extra.details=F, 
                     eps=sqrt(.Machine$double.eps))
{# The result C is the relative curvature array Bates & Watts p244 eqn. (7.16)
 # R11 is a p by p sub matrix R see Bates & Watts p236
 #  and  R2 is the m by pp sub matrix  R12 
 #                                     R22   

 p  <- dim(R11)[2]
 p2 <- dim(R11)[1]  # usually equal p except for projections
 pp <- dim(R2)[2]
 # m will usually be p+pp, but if the number of observations is small relative 
 # to pp, then Q and R2 will be truncated.
 # For smaller problems it would be quicker to do the calculation just to p+pp
 # but small problems are quick anyway.
 m  <- dim(R2)[1] 
 A  <- array(NA,c(m,p,p))  # for accel. array B&W (7.5)
 C  <- array(NA,c(m,p2,p2))
   
# R11.inv <- solve(R11)
 v <-svd(R11)
 if (all(v$d == 0)) stop("R11 is completely degenerate.")
 d <- 1/v$d
 illim <- v$d < (v$d[1] * eps)
 if (any(illim)) 
   {warning("eliminating degenerate subspace for R11.")
    d[illim] <- 0
   }
 R11.inv <- v$v %*% diag(d) %*% t(v$u) 

 R11.invT <- t(R11.inv)
 if(show.extra.details) 
      { cat("Bates and Watts matrix R11 p236)\n")
        print(R11) 
        cat("Matrix R2 (Bates and Watts  R12 \n")
        cat("                            R22 )  p236)\n")
        print(R2)
        cat("Inverse Of Matrix R11  e.g. Bates and Watts p243\n")
        print(R11.inv)
      }
  # expand R2 to make A
  # there may be a better way to do this without expanding
   k <- 0
   for(i in 1:p) 
     for(j in 1:i) 
       {k <- k + 1
        A[,i,j] <- R2[, k]
        A[,j,i] <- R2[, k]
        NULL
       }
   if(show.extra.details)
     {cat("Acceleration array (Bates and Watts matrix A p237).\n")
      cat("Parameter effects relative acceleration array:\n")
      for (i in 1:p) {print(A[i,,]); cat("\n"); NULL}
      cat("Intrinsic relative acceleration array:\n")
      for (i in (p+1):m) {print(A[i,,]); cat("\n"); NULL}
     }	

# compute C  Batts and Watts defn. (7.16) p242 and examples p244 &  p245
   for(i in 1:m)
      {C[i,,]<-(R11.invT %*% A[i,,] %*% R11.inv) * c(sqrt(p * s.sqr)); NULL}
if(show.extra.details)
  {cat("Relative curvature array (Bates and Watts matrix C p242 and eg p 243- 245).\n")
      cat("Parameter affects curvature array:\n")
      for (i in 1:p) {print(C[i,,]); cat("\n")}
      cat("Intrinsic curvature array:\n")
      for (i in (p+1):m) {print(C[i,,]); cat("\n")}
  }	

C
}


#######################################################################

#                               D matrix calculation

#######################################################################

# Notes:

#   Generating the D matrix can be computationally demanding. There are two
#   versions of the code for generating this matrix. The S version is 
#   slow and suffers from memory problems due to the way S allocates memory 
#   in loops (as of S-PLUS Version 3.0 for SUN4). The C version is
#   faster but suffers (even worse) from the memory problem.  There is a third
#   version which is written in S but using C code logic (for loops). This 
#   version is primarily for debugging the C code. Both the S and 
#   the C versions take the name of an S function 
#   as an arguement (and call it multiple times).  
#######################################################################

genD <- function(arg,...)UseMethod("genD")

genD.default <- function(func, x, d=0.01, eps=1e-4, r=6, func.args=NULL){
# This function calculates a numerical approximation of the first and second
#   derivatives (gradient and hessian) of func at the point x. The calculation
#   is done by Richardson's extrapolation (see eg. G.R.Linfield and J.E.T.Penny
#   "Microcomputers in Numerical Analysis"). The method should be used when
#   accuracy, as opposed to speed, is important.
# The result returned is a list with 3 elements as follows:
#   $D is a matrix of first(gradients) and second order partial
#      derivatives organized in the same manner as Bates and 
#      Watts. (first p columns are the gradients and the 
#      next p(p-1)/2 columns are the lower triangle of the Hessian).
#   $p is the dimension of the parameter space=dim of the tangent space.
#   $f0 is the function value at the point where the matrix D 
#        was calculated. 
#   Also the arguments (func, x, d, eps, r) are returned
#
#   modified substantially by Paul Gilbert (May, 1992)
#    from original code by Xingqiao Liu,   May, 1991.
#

#  func    name of the function. The first arg to func must be a vector.
#  x       the parameter vector.
#  func.args any additional arguments to func.
#  d       gives the fraction of x to use for the initial numerical approximation.
#  eps     is used for zero elements of x.
#  r       the number of Richardson improvement iterations.
#  v      reduction factor for Richardson iterations. This could
#      be a parameter but the way the formula is coded it is assumed to be 2.

#  Modified substantially by Paul Gilbert from first version by Xingqiao Liu.
#  This function is not optimized for S speed, but
#  is organized in the same way it is implemented in C (to facilitate checking
#  the code).

#         - THE FIRST ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F'(Xi)={F(X1,...,Xi+d,...,Xn) - F(X1,...,Xi-d,...,Xn)}/(2*d)
#                
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F''(Xi)={F(X1,...,Xi+d,...,Xn) - 2*F(X1,X2,...,Xn)+
#                     F(X1,...,Xi-d,...,Xn)}/(d^2)
#
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi, Xj IS 
#           F''(Xi,Xj)={F(X1,...,Xi+d,...,Xj+d,...,Xn)-2*F(X1,X2,...,Xn)+
#                       F(X1,...,Xi-d,...,Xj-d,...,Xn)}/(2*d^2) -
#                      {F''(Xi) + F''(Xj)}/2

#  The size of d is iteratively reduced and the Richardson algorithm is applied to
#    improve the accuracy of the computation.

   v <-2
   zt <- .001     #Rbug unix.time has trouble finding func
   # f0 <- func(x)
   f0 <- do.call(func,append(list(x), func.args))
 # zt <-unix.time(f0 <- func(x))   #  f0 is the value of the function at x.
                   # (Real value or point in sample space ie- residual vector)	
   p <- length(x)  #  number of parameters (theta)
   zt <- zt[1] *r*2*(p*(p + 3))/2
   if(zt>500) 
     cat("D matrix calculation roughly estimated to take about ",
          (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(0, length(f0),(p*(p + 3))/2)
   #length(f0) is the dim of the sample space
   #(p*(p + 3))/2 is the number of columns of matrix D.( first
   #   der. & lower triangle of Hessian)
   Daprox <- matrix(0, length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),p)
   Haprox <-  matrix(0,length(f0),r)
   for(i in 1:p)    # each parameter  - first deriv. & hessian diagonal
        {h <-h0
         for(k in 1:r)  # successively reduce h 
           {#f1 <- func(x+(i==(1:p))*h)
            #f2 <- func(x-(i==(1:p))*h) 
            f1 <- do.call(func,append(list(x+(i==(1:p))*h), func.args))
            f2 <- do.call(func,append(list(x-(i==(1:p))*h), func.args))
            Daprox[,k] <- (f1 - f2)  / (2*h[i])    # F'(i) 
            Haprox[,k] <- (f1-2*f0+f2)/ h[i]^2     # F''(i,i) hessian diagonal
            h <- h/v     # Reduced h by 1/v.
            NULL
           }
         for(m in 1:(r - 1))
            for ( k in 1:(r-m))
              {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
               Haprox[,k]<-(Haprox[,k+1]*(4^m)-Haprox[,k])/(4^m-1)
               NULL
              }
         D[,i] <- Daprox[,1]
         Hdiag[,i] <- Haprox[,1]
         NULL
        }	
   u <- p

   for(i in 1:p)   # 2nd derivative  - do lower half of hessian only
     {for(j in 1:i) 
        {u <- u + 1
            if (i==j) { D[,u] <- Hdiag[,i]; NULL}
            else 
             {h <-h0
              for(k in 1:r)  # successively reduce h 
                {#f1 <- func(x+(i==(1:p))*h + (j==(1:p))*h)
                 #f2 <- func(x-(i==(1:p))*h - (j==(1:p))*h)
                 f1 <- do.call(func, append(
                        list(x+(i==(1:p))*h + (j==(1:p))*h), func.args))
                 f2 <- do.call(func,append(
                        list(x-(i==(1:p))*h - (j==(1:p))*h), func.args))  
                 Daprox[,k]<- (f1 - 2*f0 + f2 -
                                  Hdiag[,i]*h[i]^2 - 
                                  Hdiag[,j]*h[j]^2)/(2*h[i]*h[j])  # F''(i,j)  
                 h <- h/v     # Reduced h by 1/v.
                }
              for(m in 1:(r - 1))
                 for ( k in 1:(r-m))
                   {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1); NULL}
              D[,u] <- Daprox[,1]
              NULL
             }
          }  
       }
invisible(classed(list(D=D,  # Darray constructor (genD.default)
             p=length(x),f0=f0,func=func, x=x, d=d, eps=eps, r=r), "Darray"))
}

genD.default.v2 <- function(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6){
# This function calculates a numerical approximation of the first and second
#   derivatives (gradient and hessian) of func at the point x. The calculation
#   is done by Richardson's extrapolation (see eg. G.R.Linfield and J.E.T.Penny
#   "Microcomputers in Numerical Analysis"). The method should be used 
#   accuracy, as opposed to speed, is important.
# The result returned is a list with 3 elements as follows:
#   $D is a matrix of first(gradients) and second order partial
#      derivatives organized in the same manner as Bates and 
#      Watts. (first p columns are the gradients and the 
#      next p(p-1)/2 columns are the lower triangle of the Hessian).
#   $p is the dimension of the parameter space=dim of the tangent space.
#   $f0 is the function value at the point where the matrix D 
#        was calculated. 
#
#   modified substantially by Paul Gilbert (May, 1992)
#    from original code by Xingqiao Liu,   May, 1991.
#
#  This function is not optimized for S speed, but
#  is organized in the same way it is implemented in C (to facilitate  
#  checking the code).


#  func    char string name of the function. func must have a single vector arguement.
#  x       the parameter vector.
#  d       gives the fraction of x to use for the initial numerical approximation.
#  eps     is used for zero elements of x.
#  r       the number of Richardson improvement iterations.
#  v      reduction factor for Richardson iterations. This could
#      be a parameter but the way the formula is coded it is assumed to be 2.

#         - THE FIRST ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F'(Xi)={F(X1,...,Xi+d,...,Xn) - F(X1,...,Xi-d,...,Xn)}/(2*d)
#                
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F''(Xi)={F(X1,...,Xi+d,...,Xn) - 2*F(X1,X2,...,Xn)+
#                     F(X1,...,Xi-d,...,Xn)}/(d^2)
#
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi, Xj IS 
#           F''(Xi,Xj)={F(X1,...,Xi+d,...,Xj+d,...,Xn)-2*F(X1,X2,...,Xn)+
#                       F(X1,...,Xi-d,...,Xj-d,...,Xn)}/(2*d^2) -
#                      {F''(Xi) + F''(Xj)}/2

#  The size of d is iteratively reduced and the Richardson algorithm is applied to
#    improve the accuracy of the computation.

   v <-2
   f0 <- do.call(func, append(list(x), func.args))  #  f0 is the value of the function at x.
   p <- length(x)  #  number of parameters (theta)
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(0, length(f0),(p*(p + 3))/2)
   #length(f0) is the dim of the sample space
   #(p*(p + 3))/2 is the number of columns of matrix D.( first
   #   der. & lower triangle of Hessian)
   Daprox <- matrix(0, length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),p)
   Haprox <-  matrix(0,length(f0),r)
   for(i in 1:p)    # each parameter  - first deriv. & hessian diagonal
        {h <-h0
         for(k in 1:r)  # successively reduce h 
           {f1 <- do.call(func, append(list(x+(i==(1:p))*h), func.args))
            f2 <- do.call(func, append(list(x-(i==(1:p))*h), func.args))
            Daprox[,k] <- (f1 - f2)  / (2*h[i])    # F'(i) 
            Haprox[,k] <- (f1-2*f0+f2)/ h[i]^2     # F''(i,i) hessian diagonal
            h <- h/v     # Reduced h by 1/v.
           }
         for(m in 1:(r - 1))
            for ( k in 1:(r-m))
              {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
               Haprox[,k]<-(Haprox[,k+1]*(4^m)-Haprox[,k])/(4^m-1)
              }
         D[,i] <- Daprox[,1]
         Hdiag[,i] <- Haprox[,1]
        }	
   u <- p

   for(i in 1:p)   # 2nd derivative  - do lower half of hessian only
     {for(j in 1:i) 
        {u <- u + 1
            if (i==j) { D[,u] <- Hdiag[,i] }
            else 
             {h <-h0
              for(k in 1:r)  # successively reduce h 
                {f1 <- do.call(func, append(list(x+(i==(1:p))*h + (j==(1:p))*h), func.args))
                 f2 <- do.call(func, append(list(x-(i==(1:p))*h - (j==(1:p))*h), func.args))  
                 Daprox[,k]<-(f1-2*f0+f2-Hdiag[,i]*h[i]^2-Hdiag[,j]*h[j]^2)/(2*h[i]*h[j])  # F''(i,j)  
                 h <- h/v     # Reduced h by 1/v.
                }
              for(m in 1:(r - 1))
                 for ( k in 1:(r-m))
                   Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
              D[,u] <- Daprox[,1]
             }
          }  
       }
invisible(classed(list(D=D,p=length(x),f0=f0), "Darray")) # Darray constructor (genD.default.v2)
}

genD.c <- function(func, x, d=0.01, eps=1e-4, r=6){
# See documentation in genD and genD.default 
#evalNo <<-0
#cat("evalNo = ")
   v <-2
   zt <- .001     #Rbug unix.time has trouble finding func
   f0 <- func(x)
#   zt<- unix.time(f0 <- func(x) )
   h0 <- abs(d*x)+eps*(x==0.0)
   p <- length(x)
   zt <- zt[1] *r*2*(p*(p + 3))/2
   if(zt>30) 
     cat("D matrix calculation roughly estimated to take about ",
         (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   D <- matrix(0, length(f0),(p*(p + 3))/2) 
   Daprox <-  matrix(0,length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),p)
   Haprox <-  matrix(0,length(f0),r)
   storage.mode(D)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   D<-.C("gen_D",
            list(func),             #1
            p=as.integer(length(x)),  #2
            as.double(x),
            as.double(h0),
            as.integer(length(f0)), #5
            f0=as.double(f0),       #6
            as.integer(r),
            as.integer(v),
            Haprox,        #9
            Hdiag,         #10
            Daprox,        #11
            double(length(x)),
            double(length(f0)),
            double(length(f0)),
            D=D)[c("D","p","f0")]
   invisible(classed(D, "Darray")) # Darray constructor (genD.c)
}



load.curvature.c    <- function ()
{# load C routines for use by S functions genD.c (which is NOT an improvement
 # over the S version) and R11T.c (which doesn't work).
 dyn.load(paste(DSE.HOME,"/curvature.o", sep=""))
} 

R11T.c <- function(R11.inv,p,pp){
# See documentation in R11T.s (in curvature) 
warning("usage of pp in R11T.c does not seem to be right.")
cat("this version does not work!")
   storage.mode(R11.inv)  <-"double"
   .C("R11T_c",
            R11T=matrix(0,pp,pp),
            R11.inv,
            as.integer(p),
            as.integer(pp))[["R11T"]]
}


#######################################################################

# Test routines and data for calulating curvatures in Bates and Watts.

#######################################################################



curvature.function.tests <- function( verbose=T, synopsis=T, fuzz.small=1e-14, 
	fuzz.large=1e-6, show.details=F, show.extra.details=F)
{# A short set of tests of curvature functions

  max.error <- 0
  if (synopsis & !verbose) cat("All curvature tests ...")

  if (verbose) cat("curvature test 1 ... ")


signif <- 0.05

#Rbug do.call in genD has trouble with local functions
global.assign("puromycin", function(th)
  {x <- c(0.02,0.02,0.06,0.06,0.11,0.11,0.22,0.22,0.56,0.56,1.10,1.10)
   y <- c(76,47,97,107,123,139,159,152,191,201,207,200)
   ( (th[1] * x)/(th[2] + x) ) - y
  })
global.assign("D.anal", function(th)
 {# analytic derivatives. Note numerical approximation gives a very good
  # estimate of these, but neither give D below exactly. The results are very
  # sensitive to th, so rounding error in the reported value of th could explain
  # the difference. But more likely th is correct and D has been rounded for
  # publication - and the analytic D with published th seems to work best.
  # th = c(212.70188549 ,  0.06410027) is the nls est of th for BW published D.
  x <- c(0.02,0.02,0.06,0.06,0.11,0.11,0.22,0.22,0.56,0.56,1.10,1.10)
  y <- c(76,47,97,107,123,139,159,152,191,201,207,200)
  cbind(x/(th[2]+x), -th[1]*x/(th[2]+x)^2, 0, -x/(th[2]+x)^2, 2*th[1]*x/(th[2]+x)^3)
 })
on.exit(rm(puromycin, D.anal))

# D matrix from p235. This may be useful for rough comparisons, but rounding
# used for publication introduces substanstial errors. check D.anal1-D.BW
D.BW <- t(matrix(c(
0.237812, -601.458, 0, -2.82773, 14303.4,
0.237812, -601.458, 0, -2.82773, 14303.4,
0.483481, -828.658, 0, -3.89590, 13354.7,
0.483481, -828.658, 0, -3.89590, 13354.7,
0.631821, -771.903, 0, -3.62907, 8867.4,
0.631821, -771.903, 0, -3.62907, 8867.4,
0.774375, -579.759, 0, -2.72571, 4081.4,
0.774375, -579.759, 0, -2.72571, 4081.4,
0.897292, -305.807, 0, -1.43774,  980.0,
0.897292, -305.807, 0, -1.43774,  980.0,
0.944936, -172.655, 0, -0.81173,  296.6,
0.944936, -172.655, 0, -0.81173,  296.6),   5,12))

   D.anal <- D.anal(c(212.7000, 0.0641))
#   D.anal2 <- D.anal(c(212.683, 0.0641194))
   D.calc <- genD.default("puromycin",c(212.7000, 0.0641))
#  D.calc2 <- genD.default(puromycin,c(212.683, 0.0641194))$D #est using nls
#check if col 4 is a mult. of col 2 max(abs(D.anal1-D.calc1$D))

   if (show.details)
     {cat("model A p329,data set 3 (table A1.3, p269) Bates & Watts (Puromycin example).\n")
      cat("With show.extra.details=T the calculations  on pages 235 - 252 are shown.\n")
      cat("Note Bates & Watts calculation uses s squared = 10.93^2, which is\n")
      cat(" based on both treated and untreat data. The treated data alone\n")
      cat("gives 10.39^2. This causes a slight difference in the curvature statistics.\n") 
      cat("Note that the R22 calculation using numerical gradients is quite\n")
      cat("different from the published result, but still gives similar results.\n")
      }

   r.BW   <- c(2, 12, 0.05, 0.105,        0.045,      0.21,       0.09,        1, 1.05)
   r.test <- c(2, 12, 0.05, 0.1046988133, 0.04541991, 0.21207185, 0.091999935, 1, 1.051959666)
   r.calc <-curvature.Darray(D.calc,signif=signif, 
                        show.extra.details=show.details)$stats 
   err <- r.calc - r.test
   max.error <- max(max.error, abs(err))
   ok <- all(abs(err) < fuzz.large)
   if (show.details | !ok) 
     {cat("Using numerical calculation of gradient and acceration from data:\n")
      m <- rbind(r.BW, r.test, r.calc, err)
      dimnames(m) <- list(
                    c("Bates&Watts","test value","calculated ", "difference "),
                    c(dimnames(r.calc)[[2]]))
      print(m)
      cat("Above should compared with Bates and Watts table 7.2, p 257, model A, data set 3.\n")
     }
   r.calc <- curvature.Darray(list(p=2,D=D.anal,f0=puromycin(c(212.7000, 0.0641))),
                    signif=signif, show.extra.details=show.details)$stats 
   err <- r.calc - r.test
   ok2 <- all(abs(err) < fuzz.large)
   max.error <- max(max.error, abs(err))
   if (show.details | !ok2) 
     {cat("Using D matrix of gradient and acceration calcuated with formulae Bates and Watts p234:\n")
      m <- rbind(r.BW, r.test, r.calc, err)
      dimnames(m) <- list(
                    c("Bates&Watts","test value","calculated ", "difference "),
                    c(dimnames(r.calc)[[2]]))
      print(m)
     }
   ok3 <- max(abs(D.calc$D-D.anal)) < (100 * fuzz.large)
   if (!ok3)
     {cat("max diff. between analtic and numerical D:")
      cat( max(abs(D.calc$D-D.anal)), "\n")
     }
  ok <- ok & ok2 & ok3
  all.ok <- ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("curvature test 2 ... ")

  global.assign("curvature.test.bod", function( th) 
     {x <- matrix(c(1,2,3,4,5,7),6,1)
      y <- c( 8.3, 10.3, 19.0, 16.0, 15.6, 19.8) 
      y - ( th[1] * (1 - exp( - th[2] * x)))
     })
   on.exit(rm(curvature.test.bod))
   r <- curvature.Darray(genD("curvature.test.bod",c(19.1430,  0.5311)),signif=signif, 
                show.extra.details=show.details)$stats
   r.BW   <- c(2,6, .05, 1.328,   0.184,      3.5,          0.49,        1.0, 1.01 )
   r.test <- c(2,6, .05, 1.32781, 0.18440417, 3.4990419358, 0.485941626, 1.0, 1.008372434 )
   err <- r -r.test
   max.error <- max(max.error, abs(err))
   ok <- all(abs(err) < fuzz.large )
   if (show.details | !ok) 
     {m <- rbind(r.BW, r.test, r, err)
      dimnames(m) <- list(
                    c("Bates&Watts","test value","calculated ", "difference "),
                    c(dimnames(r)[[2]]))
      print(m)
     cat("Above should compared with Bates and Watts table 7.2, p 257, model B, data set 11 (A1.4, p270.).\n")
     }
  all.ok <- all.ok & ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("curvature test 3 ... ")

  global.assign("curvature.test.isomerization", function(th)
  {x <- matrix(c(205.8,404.8,209.7,401.6,224.9,402.6,212.7,406.2,133.3,470.9,300.0,301.6
,297.3,314.0,305.7,300.1,305.4,305.2,300.1,106.6,417.2,251.0,250.3,145.1
,90.9,92.9,174.9,187.2,92.7,102.2,186.9,192.6,140.8,144.2,68.3,214.6
,142.2,146.7,142.0,143.7,141.1,141.5,83.0,209.6,83.9,294.4,148.0,291.0
,37.1,36.3,49.4,44.9,116.3,128.9,134.4,134.9,87.6,86.9,81.7,101.7
,10.5,157.1,86.0,90.2,87.4,87.0,66.4,33.0,32.9,41.5,14.7,50.2),   24,3)
   y <- c(3.541,2.397,6.6694,4.722,0.593,0.268,2.797,2.451,3.196,2.021,0.896,5.084,
5.686,1.193,2.648,3.303,3.054,3.302,1.271,11.648,2.002,9.604,7.754,11.590) 

 y - ((th[1]*th[3]*(x[,2]-x[,3]/1.632))/(1+x[,1]*th[2]+x[,2]*th[3]+x[,3]*th[4]))
})
   on.exit(rm(curvature.test.isomerization,curvature.test.bod))
   r <- curvature.Darray(genD("curvature.test.isomerization",
                           c(35.9200,0.0708,0.0377,0.1670)),signif=signif, 
                           show.extra.details=show.details)$stats
   r.BW   <-  c(4, 24, 0.05, 48.39,    0.048, 81.92, 0.081, 0.95,1.01)
#  r.test <-  c(4, 24, 0.05, 46.00941, 0.048, 81.92, 0.081, 0.95,1.01)
   r.test <-  c(4, 24, 0.05, 46.00942, 0.045942003, 77.891671, 0.077777537,
                1.0,   1.05959779)
   err <- r - r.test
   max.error <- max(max.error, abs(err))
   ok <- all(abs(err) < fuzz.large)
   if (ok) cat("Results check ok\n")
   else    cat("Results differ:\n")
   warning("curve test 3 results do not correspond with results in Bates and Watts.")
   if (show.details | !ok) 
     {m <- rbind(r.BW, r.test, r, err)
      dimnames(m) <- list(
                    c("Bates&Watts","test value","calculated ", "difference "),
                    c(dimnames(r)[[2]]))
      print(m)
#Above should compared with Bates and Watts table 7.2, p 257, model M, 
# data set 32, Bates & Watts (data table A1.5, p271).
# Above test has never worked. The intrinsic curvature array differs and 
# may need to be calculated analytically. The typical result is:
#calculated  4 24   0.05  4.600942e+01   0.045942003 
#        77.891671     0.077777537           1.00     1.05959779 
     }

  all.ok <- all.ok & ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All curvature tests completed")
     if (all.ok) cat(" OK\n")
     else    
       {cat(", some FAILED!")
        if((!is.na(max.error)) && (max.error > fuzz.small))
            cat(" ( max.error = ", max.error,")")
        cat("\n")
       }
    }
  invisible(all.ok)
}


#######################################################################

#                    end

#######################################################################


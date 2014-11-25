#eigen <- La.eigen
#svd <- function(x, nu = min(n, p), nv = min(n, p)) {
# r <- La.svd(x, nu=nu, nv=nv)
# r$v <- Conj(t(r$vt))
# r$vt <- NULL
# r
# }

version.dse <- function() cat(dse.version.information, "\n")

##############################################################################

# System dependent functions are now in the file syskern.s

# There is a short section containing installation dependent strings
# following below any general documentation which might appear next.
# (documentation may not appear here as my make procedures usually
#  filter it to a separate file.)

##############################################################################


##############################################################################

#  begin of section containing installation dependent strings.

##############################################################################


   .DSECOMPILED <- TRUE
   .DSEDUP      <- TRUE

# N.B. The default for the load.DSE.fortran function assumes the compiled object
#    is in a subdirectory named DSE.HOME/lib where DSE.HOME is a global 
#    variable set in .First.lib or in the DSE install procedure. IF THE LIBRARY 
#    IS MOVED THEN THE VALUE OF THIS VARIABLE MUST BE CHANGED.

.First.lib <- function(library, section)
  {# a previous version tried to load several functions into frame 0 to 
   # prevent problems even if the library is not attached with first=T, 
   # however, this version take the less drastic step of issuing a warning
   # if the library is not close to the beginning of the search list.

   DSE.HOME <- paste(library, "/", section, sep="")

   if (is.R())
     {if (1 < length(DSE.HOME))
         {warning(paste("package dse found in locations ", DSE.HOME,
                        "The first location is being used."))
          DSE.HOME <- DSE.HOME[1]
         }
      ok <-      require("syskern", warn.conflicts=TRUE)
      ok <- ok & require("tframe",  warn.conflicts=TRUE)
      if(!ok) warning("This package requires the syskern and tframe packages.")


      assign("DSE.HOME",DSE.HOME, envir=.GlobalEnv)            

# in very old Guide eg1.DSE.data       was called  eg1.DSE.data.diff.all.raw
# previously        eg1.DSE.data.diff  was called  eg1.DSE.data.diff.all


     } else
   if (is.S())
     {assign("DSE.HOME",DSE.HOME, where = 1)
      if (pmatch(DSE.HOME,search()) >2)
         warning("The DSE library may not work properly if it is not near the beginning of the search list. Use library(..., first=T)")
     }

   invisible( if (.DSECOMPILED) load.DSE.fortran() else TRUE )
  }



# The following variable is used in Windows dll.load
dse.win.for.tab <- c("simss" ,"smooth" ,"kfp" ,"kfprj" ,"kfepr" ,"kf" ,"simarma"
          ,"armap" ,"armaprj" ,"armaepr" ,"arma" ,"dataepr","inverse"
          ,"gend","cstat_f","efcurve_f","rlcurve_f","wep")

load.DSE.fortran <- function(from=paste(DSE.HOME,"/lib", sep=""), large=FALSE)
{if ( (0 != nchar(from)) &
     "/" != substring(from, first=nchar(from))) from <- paste(from,"/", sep="")
 if (is.R()) r <- library.dynam("dse1")  # does not use from
 else if (is.S())
    {if("windows" == Platform()$OS.type)
          {r <- dll.load(paste(from,"dsefor.dll", sep=""),dse.win.for.tab)
           if (1!=r) warning("dll load was NOT successful.")
          }
     else 
        {if (large) r <- dyn.load(paste(from,"dsefor.large.o", sep=""))
         else       r <- dyn.load(paste(from,"dsefor.o", sep=""))
    }   }
 invisible(r)
} 



help.start.DSE <- function(browser = "netscape") 
   {#browser = "netscape -mono"
    #browser = "mosaic"
    #redirecting standard input and output is necessary in order to return to S 
    if(is.R()) help.start() else if ("unix" == Platform()$OS.type())
    system(paste(browser, " ", DSE.HOME,"/help/dsehome.htm </dev/null >/dev/null &", sep=""))
    invisible()
   }

##############################################################################

#  end of section containing installation dependent strings.

##############################################################################

# Some utilities for insulating DSE code from S/R changes.

##############################################################################

if (is.R())
  { dseclass    <- class  # pre 1.4.0 .Alias(class)
   "dseclass<-" <- get("class<-")  #  .Alias(get("class<-"))
    dsescan <- function(file="",quiet=TRUE, ...){scan(file=file, quiet=quiet, ...)}
  } else
if (is.S())
  { dseclass    <- class
   "dseclass<-" <- function(x, value){ class(x) <- value ; x }
    dsescan <- function(file="",quiet=TRUE, ...){scan(file=file, ...)}
  }


##############################################################################

# The code was divided roughly into 
#   the groups listed below, but the organization has changed a little bit.
#   The actual grouping can be seen by grep ing on the string '<<<<<<' eg:
#    grep "<<<<<<" dse1.R

#    Functions which work on a model (i.e. if a model with data is allowed as
#             an arguement then the data is ignored):
#        -model summary, description, print, comparison and 
#              calculation of properties
#        -model conversion 
             
#    Functions which work on data ( and a model ):
#        -likelihood and residual calculation
#        -statistical tests
#        -parameter estimation 
#        -model reduction (this could work just on models but
#                    comparisons require data)     

#    Utility Functions:
#        -utilities for generating theoretical statistics 
#              (i.e.- model generated, not from data)
#        -data transformations
#        -model and data scaling
#        -utilities for polynomial arithmetic
#        -internal utilities used for updating objects 
#        -data interface functions

#   Of special note in the internal utilities are two programs,
#   set.arrays and set.parameters which take a model list and
#   update the representation and parameter information respectively. 
#   i.e. set.arrays uses the parameter information and ensures that the 
#    array representation is consistent while set.parameters uses the 
#    arrays and ensures that the parameter vector is consistent.

#############################################################################
#############################################################################

#functions which work on models   <<<<<<<<<<


############################################################

#     functions for model summary, description, print and
#      comparison and functions for calculation of model properties

############################################################


print.TSestModel <- function(x, ...) 
{ cat("neg. log likelihood=",x$estimates$like[1],"\n")
  if(!is.null(x$converged)) {if(!x$converged) cat(" (not converged)") }
  print(x$model,...) 
  invisible(x)
}

print.SS <- function(x, digits=options()$digits) 
    {cat("F=\n"); print(x$F,digits=digits)
     if(!is.null(x$G)) {cat("G=\n"); print(x$G,digits=digits) }
     cat("H=\n"); print(x$H,digits=digits)
     if ( is.non.innov.SS(x)) 
       {cat("Q=\n"); print(x$Q,digits=digits)
        cat("R=\n"); print(x$R,digits=digits)
        R <-x$R
       }
     if (is.innov.SS(x)) cat("K=\n"); print(x$K,digits=digits)
     if(!is.null(x$z0)) cat("initial state=\n");print(x$z0,digits=digits)
     if(!is.null(x$P0))
       cat("initial state tracking error=\n");print(x$P0,digits=digits)
     invisible(x)
    }

print.ARMA <- function(x, digits=options()$digits, L=TRUE, fuzz=1e-10) 
   # L controls the form of the display for ARMA models. 
   #If true the poly.matrix is displayed with"Ln" printed.
   {A <- x$A
     B <- x$B
     C <- x$C
     if(!is.null(x$TREND))
        cat("TREND= ", format(signif(x$TREND,digits=digits)))
     if (L)
       {cat("\nA =\n")
        for(i in 1:dim(A)[2]) 
          {for(j in 1:dim(A)[3]) 
             {cat(format(signif(A[1,i,j],digits=digits)))
              for(l in 2:dim(A)[1]) 
                 if (abs(A[l,i,j]) > fuzz)
                  {if(1==sign(A[l,i,j])) cat("+")
                   cat(format(signif(A[l,i,j],digits=digits)))
                   cat("L",l-1,sep="")
                  }
               cat("    ")
             }
           cat("\n")
         }
        cat("\nB =\n")
        for(i in 1:dim(B)[2]) 
          {for(j in 1:dim(B)[3]) 
             {cat(signif(B[1,i,j],digits=digits))
              if (2 <= dim(B)[1]) for(l in 2:dim(B)[1]) 
                 if (abs(B[l,i,j]) > fuzz)
                  {if(1==sign(B[l,i,j])) cat("+")
                   cat(signif(B[l,i,j],digits=digits))
                   cat("L",l-1,sep="")
                  }
              cat("    ")
             }
           cat("\n")
          }
        if(!is.null(x$C)) 
          {cat("\nC =\n")
           for(i in 1:dim(C)[2]) 
           {for(j in 1:dim(C)[3]) 
             {cat(signif(C[1,i,j],digits=digits))
              if (2<=dim(C)[1]) for(l in 2:dim(C)[1]) 
                 if (abs(C[l,i,j]) > fuzz)
                  {if(1==sign(C[l,i,j])) cat("+")
                   cat(signif(C[l,i,j],digits=digits))
                   cat("L",l-1,sep="")
                  }
              cat("    ")
             }
             cat("\n")
           } }
       }
     else
       {for(l in 1:dim(A)[1]) {cat("\nA(L=",l-1,")\n");print(A[l,,],digits=digits)}
        for(l in 1:dim(B)[1]) {cat("\nB(L=",l-1,")\n");print(B[l,,],digits=digits)}
        if(!is.null(x$C))
          for(l in 1:dim(C)[1]) {cat("\nC(L=",l-1,")\n");print(C[l,,],digits=digits)}
       }
     invisible(x)
} 


# summary() is a generic S function. 
# This "method" provides summary information about a model.  
 
summary.TSestModel <- function(object)
  {
   residual <- object$estimates$pred[,,drop=FALSE] - output.data(object)[,,drop=FALSE]
   sampleT <- nrow(residual)
   p <- ncol(residual)	
   Om <- t(residual) %*% residual/sampleT
   rmse <- matrix( diag(Om)^.5 ,1,p)
   dimnames(rmse) <- list(c("RMSE"), seriesNamesOutput(object))

   classed(list(  # summary.TSestModel constructor
     estimates=list(
        l=object$estimates$like[1],
        rmse=rmse,
        sampleT=sampleT,
        converged= object$converged,
        nlmin.results= (!is.null(object$nlmin.results)),
        conv.type= object$nlmin.results$conv.type,
        filter= (!is.null(object$filter)),
        smooth= (!is.null(object$smooth)) ),
     model=summary(object$model)), "summary.TSestModel")
  }


print.summary.TSestModel <- function(x, digits=options()$digits)
  {
   cat("neg. log likelihood =",x$estimates$l)
   cat("    sample length ="     ,x$estimates$sampleT, "\n")
   print(x$estimates$rmse)
   if (!is.null(x$estimates$converged))
                           cat("convergence: ",x$estimates$converged,"\n")
   if (x$estimates$nlmin.results)
                           cat("convergence type: ", x$estimates$conv.type,"\n")
   if (x$estimates$filter) cat("Includes  filter  estimates.\n")
   if (x$estimates$smooth) cat("Includes smoother estimates.\n")
   print(x$model, digits=digits)
   invisible(x)
  }



summary.SS <- function(object)
  {m <- nseriesInput(object)
   p <- nseriesOutput(object)
   n <- nrow(object$F)
   classed(list(  # summary.SS constructor
         description=object$description,
         input.series=seriesNamesInput(object),
         output.series=seriesNamesOutput(object),
         innov=is.innov.SS(object),
         m=m,
         p=p,
         n=n,
         P=n * (m+2*p),  #assumes full rank noise
         P.actual=length(coef(object)),
         constants=length(object$const),
         ICs=(!is.null(object$z0)),
         init.track=(!is.null(object$P0)) ), "summary.SS")
  }

print.summary.SS <- function(x, digits=options()$digits)
    {if (x$innov) cat("innovations form ")
     cat("state space: ")
     cat(x$description,"\n")
     cat("inputs : ", x$input.series, "\n")
     cat("outputs: ", x$output.series, "\n")
     cat("   input  dimension = ", x$m)
     cat("   state  dimension = ",x$n)
     cat("   output dimension = ",x$p,"\n")
     cat("   theoretical parameter space dimension = ",x$P,"\n")
     cat("  ",x$P.actual, " actual parameters")
     cat("   ",x$constants," non-zero constants\n")
     if (x$ICs)        cat("   Initial values specified.\n")
     else              cat("   Initial values not specified.\n")
     if (x$init.track) cat("   Initial tracking error specified.\n")
     else              cat("   Initial tracking error not specified.\n")
     invisible(x)
    }

summary.ARMA <- function(object)
  {m <- nseriesInput(object)
   p <- nseriesOutput(object)
   classed(list(  # summary.ARMA constructor
         description=object$description,
         input.series=seriesNamesInput(object),
         output.series=seriesNamesOutput(object),

         a=dim(object$A)[1]-1,
         b=dim(object$B)[1]-1,
         c=dim(object$C)[1]-1, 
         m=m,
         p=p,
         P.actual=length(coef(object)),
         constants=length(object$const),
         trend=(!is.null(object$TREND)) ), "summary.ARMA")
}
 
print.summary.ARMA <- function(x, digits=options()$digits)
    {cat("ARMA: ")
     cat(x$description,"\n")
     cat("inputs : ", x$input.series, "\n")
     cat("outputs: ", x$output.series, "\n")
     cat("     input  dimension = ", x$m)
     cat("     output dimension = ", x$p,"\n")
     cat("     order A = ", x$a)
     cat("     order B = ", x$b)
     cat("     order C = ", x$c,"\n") 
     cat("     ",x$P.actual, " actual parameters")
     cat("     ",x$constants," non-zero constants\n")
     if(x$trend) cat("     trend estimated.\n")
     else        cat("     trend not estimated.\n")
     invisible(x)
    }
 


tfplot.TSestModel <- function(x, ..., start.=NULL,end.=NULL, 
  select.inputs=NULL, select.outputs=NULL,
  Title=NULL, xlab=NULL, ylab=NULL, 
  graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
{
# plot one-step ahead estimates and actual data.
# ... is a list of models of class TSestModel.
# start. is the starting point (date) for plotting.
  model <- x
  if (is.null(Title))
     Title <- "One step ahead predictions (dotted) and actual data (solid)"
  p<- nseriesOutput(model)
  if (is.null(select.outputs)) select.outputs <-1:p
  if (all(0==select.outputs)) select.outputs <- NULL
  Ngraphs <- length(select.outputs)
  if (!is.null(select.inputs)) if (all(0==select.inputs)) select.inputs <- NULL
  m<- nseriesInput(model)
  if (is.null(m)) m <-0
  else Ngraphs <- Ngraphs+length(select.inputs)  # NULL is zero
  Ngraphs <- min(Ngraphs, graphs.per.page)
  if(reset.screen) {
     old.par <- par(mfcol = c(Ngraphs, 1), mar = mar)
     on.exit(par(old.par)) }
  names <-seriesNamesOutput(model)
  if (m!=0) names <-c(seriesNamesInput(model), names)
  if (is.null(names)) names <- rep(" ",m+p)
  if (m!=0) 
    {for (i in select.inputs) 
      {z <-NULL 
       for (model in append(list(x),list(...)) )
           z<-tbind(z,input.data(model, series=i))
       tframe(z) <-tframe(input.data(model))
       if (is.null(start.)) start.<-start(z)
       if (is.null(end.))   end.  <-end(z)
       tfOnePlot(tfwindow(z,start=start.,end=end., warn=FALSE),ylab=names[i]) # tsplot
       if(i==1) title(main = Title)
    } }
  for (i in select.outputs ) 
    {z <-c(output.data(model, series=i),
           rep(NA,dim(model$estimates$pred)[1]-periods(model$data)))
     for (model in append(list(x),list(...)))
         z<-cbind(z,model$estimates$pred[,i,drop=FALSE])
     tframe(z) <- tframe(output.data(model))
     if (is.null(start.)) start.<-start(z)
     if (is.null(end.))   end.  <-end(z)
     tfOnePlot(tfwindow(z,start=start.,end=end., warn=FALSE),ylab=names[m+i]) # tsplot
     if(i==1) title(main = Title)
    }
  invisible()
}
    


test.equal.TSestModel <- function(obj1, obj2, fuzz=0) # this could be better
  { test.equal.TSmodel( obj1$model, obj2$model, fuzz=fuzz) &
        test.equal.TSdata(obj1$data, obj2$data, fuzz=fuzz)
  }

test.equal.TSmodel <- function(obj1, obj2, fuzz=0)
{# return T if models are identical (excluding description)
  r       <- all(dseclass(obj1) == dseclass(obj2))
  if (r) r <-length(coef(obj1)) == length(coef(obj2))
  if (r) r <-all(fuzz >= abs(coef(obj1)   -  coef(obj2)))
  if (r) r <-length(obj1$location) == length(obj2$location)
  if (r) r <-all(obj1$location  ==     obj2$location)
  if (r) r <-length(obj1$i) == length(obj2$i)
  if (r) r <-all(obj1$i  ==     obj2$i)
  if (r) r <-length(obj1$j) == length(obj2$j)
  if (r) r <-all(obj1$j  ==     obj2$j)
  if (r) r <-length(obj1$const) == length(obj2$const)
  if (r) r <-all(obj1$const  ==     obj2$const)
  if (r) r <-length(obj1$const.location) == length(obj2$const.location)
  if (r) r <-all(obj1$const.location  ==     obj2$const.location)
  if (r) r <-length(obj1$const.i) == length(obj2$const.i)
  if (r) r <-all(obj1$const.i  ==     obj2$const.i)
  if (r) r <-length(obj1$const.j) == length(obj2$const.j)
  if (r) r <-all(obj1$const.j  ==     obj2$const.j)
  if (r)
    {if (is.ARMA(obj1))    r <-test.equal.ARMA(obj1,obj2, fuzz=fuzz)
     else if (is.SS(obj1)) r <-test.equal.SS(obj1,obj2, fuzz=fuzz)
    }
  r
}

test.equal.ARMA <- function(obj1, obj2, fuzz=0)
{    r <-length(obj1$l) == length(obj2$l)
     if (r) r <-all(obj1$l  ==     obj2$l)
     if (r) r <-length(obj1$const.l) == length(obj2$const.l)
     if (r) r <-all(obj1$const.l  == obj2$const.l)
     if (r) r <-length(obj1$A) == length(obj2$A)
     if (r) r <-all(fuzz >= abs(obj1$A   -     obj2$A))
     if (r) r <-length(obj1$B) == length(obj2$B)
     if (r) r <-all(fuzz >= abs(obj1$B   -     obj2$B))
     if (r) 
           {if (is.null(obj1$C)) r <- is.null(obj2$C)
             else
               {if (r) r <-length(obj1$C) == length(obj2$C)
                if (r) r <-all(fuzz >= abs(obj1$C   -     obj2$C))
           }  }
  r
}

test.equal.SS <- function(obj1,obj2, fuzz=0)
{    r <-length(obj1$F) == length(obj2$F)
     if (r) r <-all(fuzz >= abs(obj1$F   -     obj2$F))
     if (r) 
         {if(is.null(obj1$G)) r <- is.null(obj2$G)
          else
            {if (r) r <-length(obj1$G) == length(obj2$G)
             if (r) r <-all(fuzz >= abs(obj1$G   -     obj2$G))
         }  }
     if (r) r <-length(obj1$H) == length(obj2$H)
     if (r) r <-all(fuzz >= abs(obj1$H   -     obj2$H))
     if (is.innov.SS(obj1))
       {if (r) r <-length(obj1$K) == length(obj2$K)
        if (r) r <-all(fuzz >= abs(obj1$K   -     obj2$K))
       }
     else
       {if (r) r <-length(obj1$Q) == length(obj2$Q)
        if (r & (0 != length(obj2$Q)) ) r <-all(fuzz >= abs(obj1$Q - obj2$Q))
        if (r) r <-length(obj1$R) == length(obj2$R)
        if (r & (0 != length(obj2$R)) ) r <-all(fuzz >= abs(obj1$R - obj2$R))
       }
     if (r) if(is.null(obj1$z0)) r <- is.null(obj2$z0)
     else
       {if (r) r <-length(obj1$z0) == length(obj2$z0)
        if (r) r <-all(fuzz >= abs(obj1$z0   -     obj2$z0))
       }
     if (r) if(is.null(obj1$P0)) r <- is.null(obj2$P0)
     else
       {if (r) r <-length(obj1$P0) == length(obj2$P0)
        if (r) r <-all(fuzz >= abs(obj1$P0   -     obj2$P0))
       }
  r
}


McMillan.degree <- function(model,  ...) UseMethod("McMillan.degree")

McMillan.degree.TSestModel <- function(model, ...)
 {McMillan.degree(TSmodel(model),...) }

McMillan.degree.ARMA <- function(model, fuzz=1e-4, verbose=TRUE, warn=TRUE){
    z  <- roots(model, warn=warn)
    gross <- length(z)
    zz <- outer(z,z,FUN="-")
    distinct <-sum(!apply((outer(1:length(z),1:length(z),FUN="<") & (Mod(zz) <fuzz)),2,any))
    deg <- list(gross=gross,distinct=distinct)
    if (verbose)
      {cat("Assuming the model is left prime:\n")
       cat("   Without distinguishing distinct roots the degree det(A(L)) =",deg$gross,"\n")
       cat("   Distinguishing  distinct  roots       the degree det(A(L)) =",deg$distinct,"\n")
       if(!is.null(model$TREND))
          {cat("The trend adds unit roots which are added to the degree. Multiple")
           cat("unit roots are not considered distinct (but probably should be).\n")
          }
      }
    invisible(deg)
}

McMillan.degree.SS <- function(model, fuzz=1e-4){
   cat("state dimension = ",dim(model$F)[1],"/n")
   invisible()
}

stability <- function(obj, fuzz=1e-4, digits=8, verbose=TRUE) UseMethod("stability")

stability.TSestModel <- function(obj, fuzz=1e-4, digits=8, verbose=TRUE)
   {stability(TSmodel(obj), fuzz=fuzz, digits=digits, verbose=verbose)}

stability.TSmodel <- function(obj, fuzz=1e-4, digits=8, verbose=TRUE)
  {stability(roots(obj, fuzz=fuzz, randomize=FALSE),
             fuzz=fuzz, digits=digits, verbose=verbose)}

stability.roots <- function(obj, fuzz=1e-4, digits=8, verbose=TRUE) 
   {#obj <- roots(model, fuzz=fuzz, randomize=FALSE)
    s <- if (all(Mod(obj) < 1.0)) TRUE else  FALSE
    if (verbose)
      {cat("Eigenvalues of F and moduli are:\n")
       print(cbind(obj,Mod(obj)),digits=digits)
       if (s) cat("The system is stable.\n")
       else   cat("The system is NOT stable.\n")
      }
    s
   }

stability.ARMA <- function(obj, fuzz=1e-4, digits=8, verbose=TRUE) 
   {z <- roots(obj, fuzz=fuzz, randomize=FALSE)
    if (all(Mod(z) < 1.0)) s <- TRUE
    else                   s <- FALSE
    if (verbose)
      {cat("Distinct roots of det(A(L)) and moduli are:\n")
       print(cbind(1/z,Mod(1/z)),digits=digits)
       cat("\nInverse of distinct roots of det(A(L)) and moduli are:\n")
       print(cbind(z,Mod(z)),digits=digits)
       if(!is.null(obj$TREND)) cat("Trend not taken into account: ")
       if (s) cat("The system is stable.\n")
       else   cat("The system is NOT stable.\n")
      }
    s
   }


roots <- function(obj, ...) UseMethod("roots")

roots.TSestModel <- function(obj, ...){roots(TSmodel(obj),...)}

roots.SS <- function(obj, fuzz=0, randomize=FALSE) 
{   z <- eigen(obj$F)$values
    if (randomize) if (sample(c(TRUE,FALSE),1)) z <- Conj(z)
      #this prevents + - ordering of complex roots (for Monte Carlo evaluations)
    classed(z,"roots")  # constructor (roots.SS)
}


roots.ARMA <- function(obj, fuzz=0, randomize=FALSE, warn=TRUE, by.poly=FALSE) 
{   if(by.poly) z <- 1/polyroot.det(obj$A)
    else        z <- roots(to.SS(obj))
    if (fuzz!=0)
      {zz <- outer(z,z,FUN="-")  # find distinct roots within fuzz
       z <- z[ !apply((outer(1:length(z),1:length(z),FUN="<")
                 & (Mod(zz) <fuzz)),2,any)]
      }
    # add unit roots for TREND elements.
    if (!is.null(obj$TREND))
      {z <- c(rep(1,sum(0!=obj$TREND)), z)
       if (warn)
         warning("Unit roots have been added for non-zero trend elements.")
      }
    if (randomize) if (sample(c(TRUE,FALSE),1)) z <- Conj(z)
      #this prevents + - ordering of complex roots (for Monte Carlo evaluations)
    classed(z,"roots")  # constructor (roots.ARMA)
}


plot.roots <- function(x, pch='*', fuzz=0)
{if (is.TSmodel(x))    x <- roots(x, fuzz=fuzz)
 if (is.TSestModel(x)) x <- roots(x, fuzz=fuzz)
 i <- 2*pi*(1:1000)/1000
 if (max(Mod(x)) <1 )
   {plot(sin(i),cos(i),pch='.', xlab='Im', ylab='Re')
    points(Re(x),Im(x), pch=pch)
   }
 else
   {plot.default(Re(x),Im(x), pch=pch, xlab='Im', ylab='Re')
    points(sin(i),cos(i),pch='.')
   }
 points(0,0,pch='+')
 invisible(x)
}


add.plot.roots <- function(v, pch='*', fuzz=0)
{if (is.TSmodel(v))    v <- roots(v, fuzz=fuzz)
 if (is.TSestModel(v)) v <- roots(v, fuzz=fuzz)
 points(Re(v),Im(v), pch=pch)
 invisible(v)
}


observability <- function(model)  
 {# calculate singular values of observability matrix
  UseMethod("observability")
 }
observability.TSestModel <- function(model){observability(TSmodel(model)) }

observability.SS <- function(model)
{ 
FF<-    model$F
O <-    model$H
HFn <- O
for (n in 1:dim(FF)[1])  {
  HFn <- HFn %*% FF
  O <- rbind(O,HFn)
  }
svd(O)$d
}
observability.ARMA <- function(model)
{ cat("not applicable to ARMA models\n")
  invisible()
}


reachability <- function(model)  
 {# calculate singular values of reachability matrix
  UseMethod("reachability")
}

reachability.TSestModel <- function(model){reachability(TSmodel(model))}

reachability.SS <- function(model)
{FF<-    model$F
 C <-    model$G
 if (!is.null(C))
   {FnG <- C
    for (n in 1:dim(FF)[1])  
      {FnG <- FF %*% FnG
       C <- cbind(C,FnG)
      }
    cat("Singular values of reachability matrix for input: ",svd(C)$d)
   }
 if (is.innov.SS(model)) C <- model$K
 else      
  {C <- model$R
   if(dim(C)[1]==1)
     {if (any(C==0))
         {cat("State noise matrix is singular. All states are NOT excited!\n")
          return(C)
         }
      C <- 1/C
     }
   else
     {v<-svd(C)
      if (any(v$d==0))
         {cat("State noise matrix is singular. All states are NOT excited!\n")
          return(v$d)
         }
  #   C <-v$v %*% diag(1/v$d) %*% t(v$u) following is equivalent
      C <-v$v %*% (t(v$u) * 1/v$d) 
     }
   C <- model$Q %*% C
  }
 FnK <- C
 for (n in 1:dim(FF)[1])  
   {FnK <- FF %*% FnK
    C <- cbind(C,FnK)
   }
 cat("Singular values of reachability matrix for noise: ",svd(C)$d,"\n")
 invisible(svd(C)$d)
}

reachability.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


check.balance <- function(model)  
 { # calculate the difference between observability and reachability gramians 
  UseMethod("check.balance")
 }
check.balance.TSestModel <- function(model){check.balance(TSmodel(model))}

check.balance.SS <- function(model){ 
FF<-    model$F
O <-    model$H
HFn <- O
for (n in 1:dim(FF)[1])  {
  HFn <- HFn %*% FF
  O <- rbind(O,HFn)
  }
O <- t(O) %*% O # observability gramian
C <-    cbind(model$G,model$K)
FnG <- C
for (n in 1:dim(FF)[1])  {
  FnG <- FF %*% FnG
  C <- cbind(C,FnG)
  }
C <- C %*% t(C) # controllability gramian
difference <- O-C
#cat("observability gramian minus controllability gramian:\n")
#print(difference)
cat("maximum absolute difference (O-C): ", max(abs(difference)),"\n")
cat("maximum off-diagonal element of C: ", max(abs(C-diag(diag(C)))),"\n")
cat("maximum off-diagonal element of O: ", max(abs(O-diag(diag(O)))),"\n")
invisible()
}
check.balance.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


check.balanceMittnik <- function(model)  
{# calculate the difference between observability and controllability
 #   gramians with model transformed a la Mittnik 
 UseMethod("check.balanceMittnik")
}

check.balanceMittnik.TSestModel <- function(model)
   {check.balanceMittnik(TSmodel(model))}

check.balanceMittnik.SS <- function(model){ 
FF<-    model$F - model$K %*% model$H
O <-    model$H
HFn <- O
for (n in 1:dim(FF)[1])  {
  HFn <- HFn %*% FF
  O <- rbind(O,HFn)
  }
O <- t(O) %*% O # observability gramian
C <-    cbind(model$G,model$K)
FnG <- C
for (n in 1:dim(FF)[1])  {
  FnG <- FF %*% FnG
  C <- cbind(C,FnG)
  }
C <- C %*% t(C) # controllability gramian
difference <- O-C
#cat("observability gramian minus controllability gramian:\n")
#print(difference)
cat("maximum absolute difference (O-C): ", max(abs(difference)),"\n")
cat("maximum off-diagonal element of C: ", max(abs(C-diag(diag(C)))),"\n")
cat("maximum off-diagonal element of O: ", max(abs(O-diag(diag(O)))),"\n")
invisible()
}
check.balanceMittnik.ARMA <- function(model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


############################################################

#     functions for model conversion   <<<<<<<<<<

############################################################

to.SS <- function(model, ...) {UseMethod("to.SS") }
to.SS.TSestModel <- function(model, ...) 
	{l(to.SS(TSmodel(model), ...),TSdata(model))}
to.SS.SS <- function(model) {model}

to.SS.ARMA <- function(model,...)
{# convert an ARMA (or VAR) to a SS (innovations) representation
    if (is.null(model$A)) a<-0
    else a <- dim(model$A)[1] - 1  #order of polynomial arrays
    if (is.null(model$B)) b<-0
    else b <- dim(model$B)[1] - 1
    if (is.null(model$C)) cc<-0
    else cc<- dim(model$C)[1] - 1
    if ((b<=a) & (cc<=(a-1))) model <- to.SSaugment(model)
    else  model <-to.SSnested(model,...) #  (otherwise best working method) 
                  # A better approach would be an algorithm like Guidorzi's. 
 model
}


to.SSnested <- function(model, ...) {UseMethod("to.SSnested") }
to.SSnested.TSestModel <- function(model, ...)
   {to.SSnested(TSmodel(model), ...)}

to.SSnested.SS <- function(model, n=NULL, Aoki=FALSE)
{# convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <-ncol(model$F)  
  if (Aoki) return(Aoki.balance(model, n=n))
  else      return(balance.Mittnik(model, n=n)) 
}

to.SSnested.ARMA <- function(model, n=NULL, Aoki=FALSE)
{# convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <- McMillan.degree.calculation(model)$distinct
  if (Aoki) return(Aoki.balance(model, n=n))
  else      return(balance.Mittnik(model, n=n)) 
}


to.SSaugment <- function(model, ...) {UseMethod("to.SSaugment") }
to.SSaugment.TSestModel <- function(model, ...)
   {l(to.SSaugment(TSmodel(model), ...), TSdata(model))}


to.SSaugment.ARMA <- function(model, fuzz=1e-14) 
{ # convert by augmentation - state dimension may not be minimal
  # First sets A[1,,] = B[1,,] = I if that is not already the case.
   A <- model$A
   B <- model$B
   C <- model$C 
   if (fuzz  < max(abs(A[1,,]-diag(1,dim(A)[2]) )) )
      {A0.inv <- solve(A[1,,])
       A <-  polyprod(A0.inv,A)
       B <-  polyprod(A0.inv, B)
       if (!is.null(C)) C <- polyprod(A0.inv, C)
#       if (!is.null(TREND)) TREND <- A0.inv %*% TREND
       }
   if (fuzz  < max(abs(B[1,,]-diag(1,dim(B)[2]) )) )
          B<- polyprod(solve(B[1,,]), B)
   if (is.null(model$A)) a<-0
   else a <- dim(model$A)[1] - 1  #order of polynomial arrays
   if (is.null(model$B)) b<-0
   else b <- dim(model$B)[1] - 1
   if (is.null(model$C)) cc<-0
   else cc<- dim(model$C)[1] - 1
   p <- nseriesOutput(model)          # Dim of endoenous Variables.
   m <-  nseriesInput(model)          # Dim of exogenous Variables.
   if (b>a) stop("The MA order cannot exceed the AR order to convert with state augmentation.")
   if (cc>(a-1)) stop(
      " The order of the input polynomial cannot exceed the AR order -1 to convert with state augmentation.")   
  #make three parameters A,B and C have convenient order by adding 0's.
   k <- 1 + a  
#  if (b != 0)
    {BB <- array(0,c(k,dim(B)[2:3]))
     BB[1:(b+1),,] <- B
    }
   if (m!=0) 
     {CC <- array(0,c(k,dim(C)[2:3]))  
      CC[1:(cc+1),,] <- C
     }
   FF <- matrix(NA,a*p,p)
   for (i in 1:a) FF[(1+p*(i-1)):(p*i),] <- -A[a-i+2,,]
   if(a>1) FF<-cbind(rbind(matrix(0,p,(a-1)*p),diag(1,(a-1)*p)),FF)
   if (m == 0) G <-NULL
   else
     {G <- matrix(NA,a*p,m)
      for (i in 1:a) G[(1+p*(i-1)):(p*i),] <- CC[a-i+1,,] 
     }
   H <- diag(1,p)
   if(a>1) H <- cbind( matrix(0,p,(p*(a-1))),H)
   K <- matrix(NA,a*p,p)
   for (i in 1:a) K[(1+p*(i-1)):(p*i),] <- -A[a-i+2,,]+BB[a-i+2,,]
   z0 <-NULL
   if(!is.null(model$TREND))   #add a constant state which feeds into the states
     {FF<-rbind(cbind(FF,0),0) # identified with outputs (through H).
      n <-dim(FF)[1]
      FF[n,n] <-1 
      FF[n-p:1,n] <- model$TREND
      z0 <- rep(0,n)
      z0[n] <-1
      H<-cbind(H,0)
      if (m!=0) G <- rbind(G,0)
      K<- rbind(K,0)
     }                     
   descr<-c(model$description,
            " Converted to state space by state augmentation.")
   SS(F=FF,G=G,H=H,K=K,z0=z0,description=descr,
         input.names= seriesNamesInput(model),
        output.names=seriesNamesOutput(model))        
 }


gmap <- function(g, model) 
{# convert to an equivalent representation using a given matrix
 if (is.TSestModel(model)) model <- TSmodel(model)
 if  (!is.TSmodel(model)) stop("gmap expecting a TSmodel.")
 if ( is.SS(model))# transform State space model by g in GL(n)
  {n <- dim(model$F)[1]
   if (!is.matrix(g)) g<-diag(g,n) # if g is not a matrix make it into a diagonal matrix.
   if ((n !=dim(g)[1]) | (n !=dim(g)[2]) )
      stop("g must be a square matrix of dimensions equal the model state (or a scalar).")
   inv.g <- solve(g)
   model$F <-inv.g%*%model$F%*% g
   if (!is.null(model$G)) model$G <-inv.g %*%model$G
   model$H <-model$H %*% g
   if (!is.null(model$z0)) model$z0 <-c(inv.g %*%model$z0)
   if (is.innov.SS(model)) model$K <-inv.g %*% model$K
   if (is.non.innov.SS(model)) 
      {model$Q <-inv.g %*% model$Q
       model$R <-model$R
      }       
 }
 if ( is.ARMA(model))
       {# if g is not a matrix make it into a diagonal matrix.
        if (! is.matrix(g)) g<- diag(g,dim(model$A)[2]) 
	for(l in 1:dim(model$A)[1]) model$A[l,  ,  ] <- g %*% model$A[l, ,]	
	for(l in 1:dim(model$B)[1]) model$B[l,  ,  ] <- g %*% model$B[l, ,]
	for(l in 1:dim(model$C)[1]) model$C[l,  ,  ] <- g %*% model$C[l, ,]
	if(!is.null(model$TREND))   model$TREND      <- g %*%  model$TREND
       }
 set.parameters(model)
}


findg <- function(model1,model2, minf=nlmin)
{ # find the matrix which transforms between given models if it exist, 
  # otherwise the closest model...not working well
  # find g in GL(n) which minimizes the sum of squared differences between
  # parameters of models gmap(g,model1) and model2.
  # This should find the g which gives equivalence of models if that exists.
  # This procedure is rather crude and can be very slow.

  if (is.TSestModel(model1)[1]) model1 <- TSmodel(model1)
  if   (!is.TSmodel(model1)) stop("findg expecting a TSmodel.")
  if (is.TSestModel(model2)[1]) model2 <- TSmodel(model2)
  if  (!is.TSmodel(model2 )) stop("findg expecting a TSmodel.")

  if ( !is.SS(model1)| !is.SS(model2)) 
      stop("findg only works for state space models")
   n <- dim(model1$F)[1]
   if ((n!= dim(model2$F)[1])
     |(dim(model1$G)[2] != dim(model2$G)[2])
     |(dim(model1$H)[1] != dim(model2$H)[1]))
      stop("models must have the same dimensions for findg.")
   para<- c(diag(1,n))
   zzz.model1 <<- model1         # This could be done with assign(frame=1  ??)
   zzz.model2 <<- model2
   zzz.n <<-n
   func <- function(para){
      gmodel1<- gmap(matrix(para,zzz.n,zzz.n),zzz.model1)
      error <-         gmodel1$F-zzz.model2$F
      error <-c(error,(gmodel1$G-zzz.model2$G))
      error <-c(error,(gmodel1$H-zzz.model2$H))
      error <-c(error,(gmodel1$K-zzz.model2$K))
      error <-c(error,(gmodel1$R-zzz.model2$R))
      sum(error^2)
   }
   para <-minf(func,para)
   rm(zzz.model1,zzz.model2,zzz.n)
   matrix(para[[1]],n,n)
}


fix.constants <- function(model, fuzz=1e-5, constants=NULL)
{# If constants is NULL then
 # any parameters within fuzz of 0.0 or 1.0 are set to exactly 0.0 or 1.0.
 # if constants is not NULL then it should be a list with logical (T/F) arrays
 # named const.F, const.G ..., for any arrays in which there are elements
 # (not == 0 or 1) which are to be treated as constant.
 
  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("fix.constants expecting a TSmodel.")
  if (is.null(constants))
    {p <-abs(coef(model) - 1.0) < fuzz
     model$const <- c(model$const,rep(1.0,sum(p)))
     model$const.location <- c(model$const.location,model$location[p])
     model$const.i <- c(model$const.i,model$i[p])
     model$const.j <- c(model$const.j,model$j[p])
     if(is.ARMA(model)) model$const.l <- c(model$const.l,model$l[p])
     p <- (!p) & (abs(coef(model)) > fuzz) 
     model$coefficients <- coef(model)[p]
     model$location <- model$location[p]
     model$i <- model$i[p]
     model$j <- model$j[p]
     if(is.ARMA(model)) model$l <- model$l[p]
     return(set.arrays(model))
    }
  else
    return(set.parameters(TSmodel(append(model,constants))))
}


to.SS.innov <- function(model)
{ # convert to an equivalent state space innovations representation
# This assumes that the noise processes in the arbitrary SS representation are 
#   white and uncorrelated.
 if (is.TSestModel(model)) model <- TSmodel(model)
 if  (!is.TSmodel(model)) stop("to.SS.innov expecting a TSmodel.")
 if (!is.SS(model))  model <- to.SS(model)
 if ( is.non.innov.SS(model)) 
   {model$K <- model$Q %*% solve(model$R)
    model$R <- NULL
    model$Q <- NULL
   }  
 classed(model, c("innov","SS","TSmodel"))  # bypass constructor
}



to.SSOform <- function(model) {UseMethod("to.SSOform") }

to.SSOform.TSestModel <- function(model) 
  {l(to.SSOform(TSmodel(model)), TSdata(model))
  }

to.SSOform.TSmodel <- function(model)
{# convert to a SS innovations representation with a minimum number 
# of parameters by converting as much of H as possible to I matrix.
# Any remaining reductions are done by converting part of ?? to I.
# It seems there should remain n(m+2p) free parameters in F,G,H,K, and Om is 
#  determined implicitly by the residual.
 if (!is.SS(model))       model <- to.SS(model)
 if (!is.innov.SS(model)) model <- to.SS.innov(model)
 n <- dim(model$H)[2]
 p <- nseriesOutput(model)
 if (p >= n) 
   {ginv <-model$H[1:n,]
    g <-solve(ginv)
   }
 else
   {sv   <- svd(model$H)
    g    <- sv$v %*% diag(1/sv$d, ncol = length(sv$d))  %*%  t(sv$u) #right inv
    ginv <- sv$u %*% diag( sv$d,  ncol = length(sv$d))  %*%  t(sv$v) #left  inv
    g    <- cbind(g,   matrix(0,n, n-p))  # no good. these need to be full rank
    ginv <- rbind(ginv,matrix(0,n-p, n))  # and still convert H to [ I | 0 ]
    stop("This procedure is not working properly yet.") 
    # have fixed only nxp not yet nxn elements
   }
 model$F <- ginv %*% model$F %*% g
 if(!is.null(model$G)) model$G <- ginv %*% model$G
 model$H <- model$H %*% g
 model$K <- ginv %*% model$K  
 fix.constants(set.parameters(model))
}


fixF <- function(model)
{# Fix the entries in F to be constants.
 # This is a simple way to reduce the parameter dimension, but it may
 # not be a very good way to do it.
  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("fixF expecting a TSmodel.")
  if (!is.SS(model))         model <- to.SS(model)
  if ( is.non.innov.SS(model))  model <- to.SS.innov(model)
  p <-model$location == "f"
  model$const <- c(model$const, coef(model)[p])
  model$const.location <- c(model$const.location,model$location[p])
  model$const.i <- c(model$const.i,model$i[p])
  model$const.j <- c(model$const.j,model$j[p])
  p <- !p
  model$coefficients <- coef(model)[p]
  model$location <- model$location[p]
  model$i <- model$i[p]
  model$j <- model$j[p]
  cat("Remaining parameters: ",sum(p),"\n")
  if (is.null(model$G)) m <-0
  else  m <- dim(model$G)[2]
  n <- dim(model$F)[1]
  p <- dim(model$H)[1]
  cat("Theoretical parameter space dimension: ",n*(m+2*p),"\n")
  set.arrays(model)
}


to.SSChol <- function(model, ...) {UseMethod("to.SSChol") }

to.SSChol.TSestModel <- function(model, Om=NULL) 
  {if(is.null(Om)) Om <-model$estimates$cov
   l(to.SSChol(TSmodel(model), Om=Om), TSdata(model))
  }

to.SSChol.TSmodel <- function(model, Om=diag(1,nseriesOutput(model)))
{# convert to a  non.innovations SS  representation using a Cholesky 
#  decomposition of Om (the cov of the output noise). 
# Om should be an estimate of the output noise, such as returned 
#  in $estimates$cov of l.SS or l.ARMA.
# This assumes that the noise processes in the arbitrary SS representation 
# are white and uncorrelated.
 if (!is.SS(model))  model <- to.SS(model)
 if (is.innov.SS(model)) 
   {model$R <-t(chol(Om) )  # Om = RR'
    model$Q <- model$K %*% model$R
    model$K <- NULL
   }  
 classed(model, c( "non.innov","SS","TSmodel" ) )  # bypass constructor
}


to.ARMA <- function(model, ...) {UseMethod("to.ARMA") }

to.ARMA.TSestModel <- function(model, ...) 
	{l(to.ARMA(TSmodel(model), ...), TSdata(model))}

to.ARMA.ARMA <- function(model) {model}

to.ARMA.SS <- function(model, fuzz=1e-10)
{ # convert to an ARMA representation by Cayley Hamilton 
  #  (not very parsimonious)
  #ref. Aoki and Havenner, Econometric Reviews v.10,No.1, 1991, p13.
    if (is.non.innov.SS(model)) model <- to.SS.innov(model)
    FF<-model$F
    G <-model$G
    H <-model$H
    K <-model$K
    m <-dim(G)[2]
    if (is.null(m)) m <-0
    n <-dim(FF)[1]
    p <-dim(H)[1]
    poly <-  - characteristic.poly(FF) # N.B. sign change  in Aoki vs Kailath
    if ( n != length(poly))
      stop("There is some problem. The characteristic polynomial length should = state dimension.")
    A <- array(NA,c(1+n,p,p))
    A[1,,] <-diag(1,p)
    for (i in 1:n) A[i+1,,] <- diag(-poly[i],p)
    if(any(is.na(A))) stop("error in calculation of A in to.ARMA.")

#                                       i
#            Fn [i,,] corresponds to   F    
    if (n > 1)
      {Fn <- array(NA,c(n-1,n,n))
       Fn[1,,] <- FF
      }
    if (n > 2) for (i in 2:(n-1))   Fn[i,,] <- FF %*% Fn[i-1,,]

    
    HFnK <- array(NA, c(n+1, p, p))
    HFnK[1, , ] <- diag(1, p)
    HFnK[2, , ] <- H %*% K
    if (n > 1) for (i in 3:(n+1)) HFnK[i, , ] <- H %*% Fn[i-2, , ] %*% K
    B <- array(NA, c(1+n, p, p))
    B[1, , ] <- diag(1, p)
    for (i in 1:n) 
       {B[i+1, , ] <- HFnK[i+1, , ]
        for (j in 1:i) B[i+1, , ] <- B[i+1, , ] - poly[j] *  HFnK[i+1-j, , ]
       }
    if(any(is.na(B))) stop("error in calculation of B in to.ARMA.")

    if (m == 0) C <- NULL
    else
      {C <- array(NA,c(n,p,m))   
       HFnG <- array(NA,c(n,p,m)) 
       HFnG[1,,] <- H %*% G
       if (n > 1) for (i in 2:n) HFnG[i,,] <- H %*% Fn[i-1,,] %*% G
       C[1,,] <- HFnG[1,,]
       for (i in 2:n)
         {C[i,,] <- HFnG[i,,]
          for (j in 1:(i-1)) C[i,,] <- C[i,,]-poly[j]*HFnG[i-j,,]
         }
       if(any(is.na(C))) stop("error in calculation of C in to.ARMA.")
      }
 fix.constants(ARMA(A=A,B=B,C=C, input.names =  seriesNamesInput(model),
                  output.names = seriesNamesOutput(model)), fuzz=fuzz)
}


#######################################################################

#                  Utility functions  <<<<<<<<<<

############################################################

#             functions for generating  statistics  
#        from data and for generating theoretical statistics
#           (i.e.- calculated from model not data)  

############################################################

acf.M <- function(obj, ...) 
# calculate a matrix with partitions [M0|...|Mi|...|Ml] giving the cov,
#  including the exogenous series and return as first block row of Hankel.
  UseMethod("acf.M")

acf.M.TSdata <- function(obj, lag=round(6*log(periods(obj))), 
           type ="covariance", sub.mean=TRUE)
{#Estimate auto covariances from data and return as first block row of Hankel.
 # i.e. a matrix with partitions [M0|...|Mi|...|Ml] giving the cov calculated from 
 #  the data, including the exogenous series and return as first block row of Hankel.
 #  Each Mi is a p by (m+p) matrix, (m is the dimension of the exogeneous 
 #  series and p is the dimension of endogeneous series)
 #  ie.  Mi = [ cov{y(t),y(t-i)} | cov{y(t),u(t-i)} ] 
 #  N.B. - The part of the first block corresponding to y(t)y(t) may need to be discarded
 #     for Aoki's technique.  This will reverse the order of y and u!!!!
 # if type == "correlation" the result is scaled to give autocorrelations.
  data <- freeze(obj)
  p <- ncol(output.data(data))
  m <- nseriesInput(data)
  if (is.null(m)) m <-0
  d <- cbind(output.data(data),input.data(data))
  sampleT <-periods(data)
  if (sub.mean) d <- d- t(matrix(apply(d,2,mean),dim(d)[2],sampleT))
#  if (type == "correlation") d <- d %*% diag(1/diag(t(d)%*%d/sampleT)^.5)
#  z <-acf(d, type = type, plot=FALSE)$acf
  z <- array(NA, c(lag,p,dim(d)[2]))
  for (i in 1:lag) 
    {z[i,,] <- (t(d[i:sampleT,1:p]) %*% d[1:(sampleT+1-i),]) / (sampleT+1-i)
    }
  M <- NULL
  for (i in 1:dim(z)[1])  M <- cbind(M, z[i,1:p,])
  if (type == "correlation")
    {ro <- matrix(1/diag(M)^.5,p,p)         # this will mess up if m != 0
     if (m!=0) cat("input variables not yet handled correctly!")
     M  <- array(ro*t(ro),dim(M)) * M
    }
  M
}

acf.M.TSestModel <- function(obj, ...) {acf.M(TSmodel(obj), ...)}

acf.M.TSmodel <- function(obj, lag=NULL, type ="covariance", Psi=NULL)
{# Construct a matrix with partitions [M0|...|Mi] giving the theoretical cov calculated from 
 #  the model, including the exogenous parameters and return as first block row of Hankel.
 #  Each Mi is a p by (m+p) matrix, (m is the dimension of the exogeneous 
 #  series and p is the dimension of endogeneous series)
 #  ie.  Mi = [ cov{y(t),y(t-i)} | cov{y(t),u(t-i)} ] 
 #  The innovations cov Psi could depend on the model provided but does not yet.
 #  If specified it is used. If not specified it is set to I.
 #  N.B. - The part of the first block corresponding to Gamma0= E{y(t)y(t)'} 
 #    The first block of the result (Gamma0) may need to be discarded for Vaccaro's and
 #     for Aoki's technique.  This will reverse the order of y and u!!!!
 # if type == "correlation" the result is scaled to give autocorrelations.
 
if (!is.SS(model)) model <- to.SS(obj) else model <- obj
if(is.null(Psi)) Psi <- diag(1,dim(model$H)[1])
if(is.null(lag)) lag <- 3*dim(model$F)[1]
cat (" Warning: Cov generation has not been tested(and doesn't work).\n")
if ( is.ARMA(model))
  #  M ={ Mi }={ Ci-1|Bi| -Ai}, i=2,...,k. k=max(a,b,cc). Assumes I=A[1,,]=B[1,,]
  {A <- model$A
   B <- model$B
   C <- model$C
   a <- dim(A)[1] - 1      # order of polynomial arrays
   if (is.na(a)) a <- 0
   b <- dim(B)[1] - 1
   if (is.na(b)) b <- 0
   cc <- dim(C)[1] - 1 
   if (is.na(cc)) cc <- 0
   m <- dim(C)[3]          # Dim of exogenous Variables.
   if (is.null(m))  m <- 0                         
   #make three parameters A,B and C have convenient order by adding 0's.    
   k <- 1 + max(a,cc,b)  
   if (a != 0) 
    {AA <- array(0,c(k,dim(A)[2:3]))
     AA[1:(a+1),,] <- A
    }
   if (b != 0)
    {BB <- array(0,c(k,dim(B)[2:3]))
     BB[1:(b+1),,] <- B
    }
   if (m != 0)
    {CC <- array(0,c(k,dim(C)[2:3]))  
     CC[1:(cc+1),,] <- C
    }
   M <- NULL
cat (" Warning: Cov generation does not yet work for ARMA models.\n")
   for(i in 2:k)  
       {if (m != 0) M <- cbind(M,CC[(i-1),,]) 
        if (b != 0) M <- cbind(M, BB[i,,])      # constant term ignored (=I)
        if (a != 0) M <- cbind(M,-AA[i,,])      # constant term ignored (=I)
       }
   if(!is.null(model$TREND))
      cat(" Warning: Theoretical cov generation does not account for trends.")
  }   
if ( is.SS(model))
   {FF<-model$F
    G <-model$G
    H <-model$H
    if (is.innov.SS(model)) K <-model$K
    else          K <- model$Q %*% solve(model$R)
    P  <- Riccati(FF,K %*% Psi %*%t(K))  # V&V (6)
    M <- H %*% P %*% t(H) + Psi  # Gamma0      V&V (7)
    Om <-  K %*% Psi + FF %*% P %*% t(H) # V&V (7)
    FnKG <- FF %*% P %*% t(H) + K %*% Psi    # V&V (5) 
    if (!is.null(G)) 
      {FnKG <- cbind(FnKG, G)       # + G is NOT correct
       M <- cbind(M,G)
      }
    i<-0    # M should have at least 2 blocks or Hankel shift does not work.
    while(i<=lag)
       {FnKG <- FF %*% FnKG
        M <- cbind(M,H %*% FnKG)    
        i <-i+1   
       }
  }    
if (type == "correlation")
    {p  <- dim(H)[1]
     if (!is.null(G)) cat("input variables not yet handled correctly!")
     ro <- matrix(1/diag(M)^.5,p,p)                  # this will mess up if m != 0
     M  <- array(ro*t(ro),dim(M)) * M
    }
list(M=M, Om=Om,P=P)       
}



Riccati <- function(A, B, fuzz=1e-10, iterative=FALSE)
{warning("This procedure has not been tested!")
 if (!iterative) 
  {n <- dim(A)[1]
   Atinv <-solve(t(A))    # symplectic matrix Vaughan (10)(12), R=0
   S   <- rbind(cbind(Atinv     , diag(0,n)),
               cbind(B %*% Atinv, A)) 
   Q <- eigen(S)
   Q <- Q$vectors[,rev(order(Mod(Q$values)))]   # This may have imaginary parts.
   P <- Re( Q[(n+1):(n+n),1:n] %*% solve(Q[1:n,1:n]) ) #This should not have any significant im parts.
  }
 else
  {P<- diag(0,dim(A)[1])
   i <-0
   repeat    # yuk
     {P <- A%*% P %*% t(A) + B
      i <- i+1
      if (i>1000) break
      if (fuzz > max(abs(P-A %*% P %*% t(A) - B))) break
  }  }
 if (fuzz < max(abs(P-A %*% P %*% t(A) - B)))
      warning("Riccati failed! Result does not solve the Riccati equation!")
 P
}


markov.parms <- function(model, blocks=NULL) 
{if (is.TSestModel(model)) model <- TSmodel(model)
 if  (!is.TSmodel(model)) stop("markov.parms expecting a TSmodel.")
 if ( is.ARMA(model))
  #  M ={ Mi }={ Ci-1|Bi| -Ai}, i=2,...,k. k=max(a,b,cc). Assumes I=A[1,,]=B[1,,]
  {A <- model$A
   B <- model$B
   C <- model$C
   a <- dim(A)[1] - 1      # order of polynomial arrays
   if (is.na(a)) a <- 0
   b <- dim(B)[1] - 1
   if (is.na(b)) b <- 0
   cc <- dim(C)[1] - 1 
   if (is.na(cc)) cc <- 0
   m <- dim(C)[3]          # Dim of exogenous Variables.
   if (is.null(m))  m <- 0                         
   #make three parameters A,B and C have convenient order by adding 0's.    
   if (is.null(blocks))
      blocks <- 1 + max(2,a,cc,b) 
    # if blocks is not at least 3 the Hankel shift does not work
   if (a != 0) 
    {AA <- array(0,c(blocks,dim(A)[2:3]))
     AA[1:(a+1),,] <- A
    }
   if (b != 0)
    {BB <- array(0,c(blocks,dim(B)[2:3]))
     BB[1:(b+1),,] <- B
    }
   if (m != 0)
    {CC <- array(0,c(blocks,dim(C)[2:3]))  
     CC[1:(cc+1),,] <- C
    }
   M <- NULL
   if (b != 0) cat (" Warning: This has only been developed for SS and VARX models.\n")
   if(!is.null(model$TREND)) cat(" Warning: Markov parameter generation does not account for trends.")
   for(i in 2:blocks)  
       {if (m != 0) M <- cbind(M,CC[(i-1),,]) 
        if (b != 0) M <- cbind(M, BB[i,,])      # constant term ignored (=I)
        if (a != 0) M <- cbind(M,-AA[i,,])      # constant term ignored (=I)
       }
  }   
else if ( is.SS(model))
   {FF<-model$F
    G <-model$G
    H <-model$H
    if (is.innov.SS(model)) K <-model$K
    else               K <- model$Q %*% solve(model$R)
    if (is.null(blocks)) blocks <- 1+dim(FF)[1]
    FF <- FF - K %*% H  # model transformed a la Mittnik so lagged outputs are inputs
    FnGK <- cbind(G,K)
    M <- H %*% FnGK
    i<-0 # M should have at least 2 blocks or Hankel shift does not work.
    stopp <- FALSE
    while((i<=blocks) & !stopp)   #  no. of blocks affect Hankel size
       {FnGK <- FF %*% FnGK
        M <- cbind(M,H %*% FnGK)
        i <-i+1   # count should not be necessary, but insures an end.
        stopp <- (i>3) & ( max(abs(FnGK)) < 1e-15)
       }
  }
else stop("markov.parms requires an ARMA or SS model.")
M       
}



############################################################

#     polynomial utility functions   <<<<<<<<<<

############################################################


polyprod <- function(a,b)
{ # product of two polynomials.
# The convention used is by poly.value and polyroot is constant first,
#  highest order coef. last. The reverse convention could also be used for multiplication.
# This function handles scalar (ie. non-matrix) and matrix polynomials.
# Scalar polynomials are vectors of length 1+the polynomial order.
# Polynomial matrices are defined as 3 dimensional arrays with the last 2
# dimensions as the matrix dimension and the first equal 1+the
# polynomial order.

    pprod <- function(a,b)  # local function, product of non-matrix polys.
       {n <- length(a) +length(b) -1
        if (is.null(a))    return(NA)
        if (is.null(b))    return(NA)
        if (0 ==length(a)) return(NA)
        if (0 ==length(b)) return(NA)
        if (any(is.na(a))) return(NA)
        if (any(is.na(b))) return(NA)
	r <- rep(NA, n)
        z <- outer(a, b) 
        zi <- 1 + outer(1:length(a)-1,1:length(b)-1,"+")
	for(i in 1:n) r[i]<- sum(z[i==zi])
	r
       }
   psum <- function(a,b)  # local function, sum of non-matrix polys.
    {if (length(a) < length(b)) return(c(a,rep(0,length(b)-length(a))) + b)
     else                       return(c(b,rep(0,length(a)-length(b))) + a)
    }
   if (is.vector(b)  && (is.array(a) | is.matrix(a)))
     {z <- b; b <- a; a <- z } # scalar multiplication commutes (even for scalar polynomials)
   if (is.null(a))      r <- NULL  
   else if (is.null(b)) r <- NULL   
   else if (is.vector(a))  
          {if      (is.vector(b)) r <-pprod(a,b)
           else if (is.matrix(b))
              {r <- array(NA,c(length(a),dim(b)))
               for (i in 1:(dim(b)[1])) 
                  for (j in 1:(dim(b)[2]))
                     r[,i,j] <- pprod(a,b[i,j])
              }
           else if (is.array(b))
              {r <- array(NA,c(length(a)+dim(b)[1]-1,dim(b)[2:3]))
               for (i in 1:(dim(b)[2])) 
                  for (j in 1:(dim(b)[3]))
                     r[,i,j] <- pprod(b[,i,j],a)
              }
          }
   else if (is.matrix(a))
          {if (is.matrix(b)) r <- a %*% b
           else if (is.array(b))
              {if (dim(a)[2] != dim(b)[2]) 
                  stop("Matrix polynomial dimensions do not conform.")
               r <- array(0,c(dim(b)[1],dim(a)[1], dim(b)[3]))
               #for (i in 1:(dim(a)[1])) for (j in 1:(dim(b)[3]))
                 #for (k in 1:(dim(a)[2])) r[,i,j] <- psum(r[,i,j], pprod(a[i,k],b[,k,j]))
               for (k in 1:(dim(b)[1]))
                 r[k,,] <- a %*% array(b[k,,],dim(b)[2:3])
                 #array above insures b a matrix, drop=FALSE cannot be used
              }
          }
   else if (is.array(a))
        if (is.matrix(b))
          {if (dim(a)[3] != dim(b)[1])
              stop("Matrix polynomial dimensions do not conform.")
           r <- array(0,c(dim(a)[1],dim(a)[2], dim(b)[2]))
           #for (i in 1:(dim(a)[2])) for (j in 1:(dim(b)[2]))
           #  for (k in 1:(dim(a)[3]))  r[,i,j] <- psum(r[,i,j], pprod(a[,i,k],b[k,j]))
           for (k in 1:(dim(a)[1]))
                 r[k,,] <- array(a[k,,],dim(a)[2:3]) %*% b
                 #array above insures b a matrix, drop=FALSE cannot be used
          }
        else if (is.array(b))
          {if (dim(a)[3] != dim(b)[2]) 
              stop("Matrix polynomial dimensions do not conform.")
           r <- array(0,c(dim(b)[1]+dim(a)[1]-1,dim(a)[2], dim(b)[3]))
           for (i in 1:(dim(a)[2])) 
              for (j in 1:(dim(a)[3]))
                 for (k in 1:(dim(a)[3]))
                   r[,i,j] <- psum(r[,i,j], pprod(a[,i,k],b[,k,j]))
          }
   else stop("polynomial product not defined for these objects")
r
}



polysum <- function(a,b){
   # sum of two polynomials (including polynomial arrays)
    psum <- function(a,b)  # local function for non-matrix polys.
 {if (length(a) < length(b))  {r <- b;  r[1:length(a)] <-r[1:length(a)] + a }
  else  {r <- a; r[1:length(b)] <-r[1:length(b)] + b }
  r 
 }
   if (is.null(a) )     r <- b
   else if (is.null(b)) r <- NULL
   else if (is.vector(a) && is.vector(b)) r <-psum(a,b)
   else if (is.array(a) )
        if (is.array(b))
           if ( all(dim(a)[2:3] == dim(b)[2:3]))
             {if (dim(a)[1] < dim(b)[1])
                    {r <- b
                     r[1:(dim(a)[1]),,] <-r[1:(dim(a)[1]),,] + a
                    }
                 else
                    {r <-  a
                     r[1:(dim(b)[1]),,] <-r[1:(dim(b)[1]),,] + b
                    }
             }
           else stop("polynomial matrix dimensions must agree")
        else if (is.vector(b))
          {r <- array(NA,c(max(length(b),dim(a)[1]),dim(a)[2:3]))
           for (i in 1:(dim(a)[2])) 
              for (j in 1:(dim(a)[3]))
                 r[,i,j] <- psum(a[,i,j],b)
          }
   else if (is.vector(a) && is.array(b))
          {r <- array(NA,c(max(length(a),dim(b)[1]),dim(b)[2:3]))
           for (i in 1:(dim(b)[2])) 
              for (j in 1:(dim(b)[3]))
                 r[,i,j] <- psum(b[,i,j],a)
          }
   else stop("polynomial sum not defined for these objects")
r             
}

polyroot.det <- function(A)
{#roots of the determinant of A.  Note: polydet is slow. There is room for improvement here!
z <- polydet(A)
if (length(z)==1) 
  stop(paste("root cannot be calculated for degree 0 determinant = ", as.character(z)))
polyroot(z)
}

polydet <- function(A)
{# recursive pessimist: Life is just one damned thing before another.
 # attributed to Richard Bird in The Mathematical Intelligencer, Winter 1994.
 #Recursively form the determinant of a polynomial matrix A, where the first
 #  dim stores the coefs. of the polynomial for each A[,i,j].
	n <- dim(A)[2]  
        if (n != dim(A)[3])
           stop( "The determinant is only defined for square polynomial matrices")
	if(1 == n) r <- c(A)
	else 
          {r<- 0   # previously NULL
           for (i in 1:n) 
             {if(!all(0==A[,i,1]))#not nec.but faster for sparse arrays
                   r<- polysum(r,(-1)^(i+1)*
          polyprod(A[,i,1],Recall(A[,(1:n)[i!=(1:n)],2:n,drop=FALSE])))
              r[is.na(r)] <- 0
              if (any(0==r)) 
                   {if (all(r==0)) r <- 0
                    else  r <-r[0 == rev(cumprod(rev(r==0)))] #remove trailing zeros
                    if(0==length(r)) r <- 0
             }     } 
          }
r
}

poly.value <- function(coef, z)
{# evaluate a polynomial given by coef (constant first) at z
 # could be extended for matrix coef.
#           n-1           n-2
#  coef[n]*z   + coef[n-1]*z   + ... + coef[1]  
  coef %*% z^(0:(length(coef)-1))
}


characteristic.poly <- function(A){
	# coefficients of the characteristic polynomial of a matrix
	# return a vector of the coefficients of the characteristic polynomial of A.
	# ref. Kailath "Linear Systems" p657

   tr <- function(A){ # calc the trace of a matrix
     sum(diag(A))
   }
   n <- dim(A)[1]
   if (n != dim(A)[2]) stop(" arguement must be a square matrix.")
   s <-array(0,c(n,n,n))
   if (n==1) a <- -A
   else
     {a <- rep(0,n)  
      s[1,,] <- diag(1,n)
      for (i in 1:(n-1))
        {a[i] <- -tr(s[i,,]%*%A)/i
         s[i+1,,] <- (s[i,,]%*%A) + diag(a[i],n)
        }
      a[n] <- -tr(s[n,,]%*%A)/n
     }
   a
}
companion.matrix <- function(A)
{# return the (top) companion matrix for a 3 dim array A (polynomial matrix), 
#  where the 1st dim corresponds to coefs. of polynomial powers (lags).
# ref. Kailath "Linear Systems" p659
  p <- dim(A)[2]
  if (p!= dim(A)[3]) 
    stop("companion matrix can only be computed for square matrix polynomials")
  l <- dim(A)[1]  # 1+ order of A
  C <- rbind(matrix(0,p,l*p), cbind(diag(1,(l-1)*p,(l-1)*p), matrix(0,(l-1)*p,p)))
  for (i in 1:l) C[1:p,((l-1)*p+1):(l*p)] <- -A[l,,]
  C
}

############################################################

#     internal utility functions    <<<<<<<<<<

############################################################


read.int <- function(prmt)
   {err <-TRUE
    while (err)
     {cat(prmt)
      n <- as.integer(readline()) # crude. this truncates reals
      if (is.na(n)) cat("value must be an integer\n")
      else err <- FALSE
     }
    n
   }


nseriesInput <- function(x)UseMethod( "nseriesInput")
nseriesInput.default <- function(x)
   if (is.null(x$input)) 0 else nseries(x$input)
nseriesInput.SS <- function(x)   {if (is.null(x$G)) 0 else  dim(x$G)[2] }
nseriesInput.ARMA <- function(x) {if (is.null(x$C)) 0 else dim(x$C)[3] }
nseriesInput.TSestModel <- function(x){nseriesInput(x$data)}

nseriesOutput <- function(x)UseMethod("nseriesOutput")
nseriesOutput.default <- function(x)
   if (is.null(x$output)) 0 else nseries(x$output)
nseriesOutput.SS <- function(x){dim(x$H)[1] }
nseriesOutput.ARMA <- function(x){dim(x$A)[2] }
nseriesOutput.TSestModel <- function(x){nseriesOutput(x$data)}


check.consistent.dimensions <- function(obj, ...)
   { # check data & model dimensions
    UseMethod("check.consistent.dimensions")
   }
check.consistent.dimensions.TSdata <- function(obj, model)
   {check.consistent.dimensions(model, obj)}
check.consistent.dimensions.TSestModel <- function(obj, data)
   {if(missing(data)) data <- obj$data
    check.consistent.dimensions(obj$model,data)
   }

check.consistent.dimensions.SS <- function(obj,data=NULL)
 {m <-dim(obj$G)[2]
  n <-dim(obj$F)[1]
  p <-dim(obj$H)[1]
  if (n!= dim(obj$F)[2]) stop("Model F matrix must be square.")
  if (n!= dim(obj$H)[2])
      stop("Model H matrix have second dimension consistent with matrix F.")
  if (!is.null(obj$G)) if(n!= dim(obj$G)[1])
      stop("Model G matrix have first dimension consistent with matrix F.")
  if (!is.null(obj$K)) if(n!= dim(obj$K)[1])
      stop("Model K matrix have first dimension consistent with matrix F.")
  if (!is.null(obj$K)) if(p!= dim(obj$K)[2])
      stop("Model K matrix have second dimension consistent with matrix H.")

  if (!is.null(data))
   {if(dim(obj$H)[1] != nseriesOutput(data))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(obj$G))
      {if(0 != nseriesInput(data))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(obj$G)[2] != nseriesInput(data))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(TRUE)
 }

check.consistent.dimensions.ARMA <- function(obj,data=NULL)
 {p <-dim(obj$A)[2]
  if (p!= dim(obj$A)[3]) stop("Model A array dim 2 and 3 should be equal.")
  if (p!= dim(obj$B)[2]) stop("Model B array dim inconsistent with array A.")
  if (p!= dim(obj$B)[3]) stop("Model B array dim inconsistent with array A.")
  if (!is.null(obj$C))
      if (p!= dim(obj$C)[2]) stop("Model C array dim inconsistent with array A.")

  if (!is.null(data))
   {if(dim(obj$A)[2] != nseriesOutput(data))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(obj$C))
      {if(0 != nseriesInput(data))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(obj$C)[3] != nseriesInput(data))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(TRUE)
 }

check.consistent.dimensions.default <- function(obj, ...)
   {stop(paste("No method for object of class ", dseclass(obj), "\n"))}



TSestModel <- function(obj) UseMethod("TSestModel")
TSestModel.TSestModel <- function(obj) {obj} # return everything


TSmodel <- function(obj, ...) UseMethod("TSmodel")
TSmodel.TSmodel <- function(obj){obj}

TSmodel.TSestModel <- function(obj)
  {# Return a TSmodel object but also retains convergence info (if not null).
   obj$model$converged <- obj$converged
   obj$model
  }




ARMA <- function(A=NULL, B=NULL, C=NULL, TREND=NULL, description=NULL,
          names=NULL, input.names=NULL, output.names=NULL) 
  {if  (is.null(A)) stop("specified structure is not correct for ARMA model.")
   # and fix some simple potential problems
   if(is.null(dim(A)))   A <- array(A, c(length(A),1,1))
   if(is.null(B)) stop("B array must be specified for ARMA class models.")
   if(is.null(dim(B)))   B <- array(B, c(length(B),1,1))
   if(2==length(dim(B))) B <- array(B, c(1, dim(B)))
   model <- list(A=A, B=B, C=C, TREND=TREND, description=description)
   dseclass(model) <- c("ARMA","TSmodel") # constructor
   if(!is.null(names)) seriesNames(model) <- names
   else
     {if(!is.null( input.names))  seriesNamesInput(model) <- input.names
      if(!is.null(output.names)) seriesNamesOutput(model) <- output.names
     }
   check.consistent.dimensions(model, data=NULL)
   set.parameters(model)
  }



SS <- function(F.=NULL, G=NULL, H=NULL, K=NULL, Q=NULL, R=NULL, z0=NULL, P0=NULL,
             description=NULL,
	     names=NULL, input.names=NULL, output.names=NULL)   
  {if (is.null(F.) | is.null(H))
       stop("specified stucture is not correct for SS model.")
   # and fix some simple potential problems
   if(1==length(F.))  
     {if(!is.matrix(F.))                 F. <- matrix(F.,1,1)
      if(!is.matrix(H))                  H  <- matrix(H,length(H),1)
      if(!is.null(G) && !is.matrix(G))   G  <- matrix(G,1,length(G))
      if(!is.null(K) && !is.matrix(K))   K  <- matrix(K,1,length(K))
      if(!is.null(Q) && !is.matrix(Q))   Q  <- matrix(Q,1,length(Q))
     }
   model <- list(F=F., G=G, H=H, K=K, Q=Q, R=R, z0=z0, P0=P0,
                 description=description)
   if (!is.null(model$K)) dseclass(model) <- c("innov","SS","TSmodel" ) else
   if (!is.null(model$Q)) dseclass(model) <- c( "non.innov","SS","TSmodel") # constructor
   else stop("specified stucture is not correct for SS model.")
   if(!is.null(names)) seriesNames(model) <- names
   else
     {if(!is.null( input.names))  seriesNamesInput(model) <-  input.names
      if(!is.null(output.names)) seriesNamesOutput(model) <- output.names
     }
   check.consistent.dimensions(model, data=NULL)
   set.parameters(model)
  }



"coef<-" <- function(object, value) UseMethod("coef<-")
"coef<-.default" <- function(object, value) {
   object$coefficients <- value
   object
   }
   
coef.TSmodel <- function(object) object$coefficients
coef.TSestModel <- function(object) coef(TSmodel(object))

is.TSmodel <- function(obj){inherits(obj,"TSmodel")}
is.TSestModel <- function(obj){inherits(obj,"TSestModel")}
is.SS <- function(obj){inherits(obj,"SS")}
is.innov.SS <- function(obj){inherits(obj,"SS")& inherits(obj,"innov")}
is.non.innov.SS <- function(obj){inherits(obj,"SS")&inherits(obj,"non.innov")}
is.ARMA <- function(obj){inherits(obj,"ARMA")}

set.parameters <- function(model)  
  { # complete parameter info. based on representation info. 
   set.parameters.TSmodel(TSmodel(model))
  }

set.parameters.TSmodel <- function(model)  
   UseMethod("set.parameters.TSmodel")

set.parameters.TSmodel.SS <- function(model) { 

 locateSS <- function(A,Ac,a,I,J,plist)# local function for locating parameters
  {indicate <-  (A==1.0)                # constants
   if (!is.null(Ac)) indicate <- indicate | Ac  # Ac is T for fixed entries
   plist$const <- c(plist$const,A[indicate])
   plist$const.location <- c(plist$const.location,rep(a,sum(indicate)))
   if (a!="z")
     {plist$const.i <- c(plist$const.i,row(A)[indicate]) 
      plist$const.j <- c(plist$const.j,col(A)[indicate])
     }
   else
     {plist$const.i <- c(plist$const.i,(1:length(A))[indicate]) 
      plist$const.j <- c(plist$const.j,(1:length(A))[indicate]) #dup.but ok
     }
   indicate <- (A!=0.0) & (A!=1.0) 
   if (!is.null(Ac)) indicate <- indicate & (!Ac)
   coef(plist) <- c(coef(plist), A[indicate])          # parameters
   plist$location <- c(plist$location,rep(a,sum(indicate)))
   if (a!="z")
     {plist$i <- c(plist$i,row(A)[indicate])
      plist$j <- c(plist$j,col(A)[indicate])
     }
   else
     {plist$i <- c(plist$i,(1:length(A))[indicate]) 
      plist$j <- c(plist$j,(1:length(A))[indicate]) #dup.but ok
     }
   plist
 }# end locateSS    

    m <-dim(model$G)[2]
    n <-dim(model$F)[1]
    if (n!= dim(model$F)[2]) stop("Model F matrix must be square.")
    if (n!= dim(model$H)[2])
      stop("Model H matrix have second dimension consistent with matrix F.")
    if (!is.null(model$G)) if(n!= dim(model$G)[1])
      stop("Model G matrix have first dimension consistent with matrix F.")
    p <-dim(model$H)[1]
    plist <- locateSS(model$F, model$const.F,"f",n,n,
                 list(coef=NULL,location=NULL,i=NULL,j=NULL,
                 const=NULL,const.location=NULL,const.i=NULL,const.j=NULL))
    if(!is.null(m)) plist <- locateSS(model$G,model$const.G,"G",n,m,plist)
    plist <- locateSS(model$H,model$const.H,"H",p,n,plist)
    if(!is.null(model$z0)) plist <- locateSS(model$z0,model$const.z0,"z",p,n,plist)
    if(!is.null(model$P0)) plist <- locateSS(model$P0,model$const.P0,"P",p,n,plist)
    if (is.innov.SS(model)) 
      {plist <- locateSS(model$K,model$const.K,"K",n,p,plist)
       # note const.H, etc are logical arrays (if not NULL) to indicate
       #      parameters which are to remain fixed, so that set.parameters
       #       knows to put them in const. This feature has not been used much.
      }
    else
      {plist <- locateSS(model$Q,model$const.Q,"Q",n,n,plist)
       plist <- locateSS(model$R,model$const.R,"R",p,p,plist)
      }

#    list.add(model, names(plist) ) <- plist
    model[ names(plist) ] <- plist
    model
} #end set.parameters.SS

set.parameters.TSmodel.ARMA <- function  (model) { 

 locateARMA <- function(A,a,I,J,L,plist){ # local function for locating parameters
  ind <- function(x, i) # equivalent of row and col for higher dim arrays.
   {
	d <- dim(x)
	id <- 1:length(d)
	id <- c(i, id[i!=id])
	y <- array(1:(d[i]), dim(x)[id])
	aperm(y, order(id))
   }
   indicate <-  (A==1.0)                # constants
   plist$const <- c(plist$const,A[indicate])
   plist$const.location <- c(plist$const.location,rep(a,sum(indicate)))
   if (a!="t")
     {plist$const.l <- c(plist$const.l,ind(A,1)[indicate]) 
      plist$const.i <- c(plist$const.i,ind(A,2)[indicate]) 
      plist$const.j <- c(plist$const.j,ind(A,3)[indicate])
     }
   else    # trend is a vector not an array 
     {plist$const.l <- c(plist$const.l,rep(0,sum(indicate))) 
      plist$const.i <- c(plist$const.i, (1:length(A))[indicate]) 
      plist$const.j <- c(plist$const.j,rep(0,sum(indicate)))
     }
   indicate <- (A!=0.0) & (A!=1.0)
   coef(plist) <- c(coef(plist), A[indicate])          # parameters
   plist$location <- c(plist$location,rep(a,sum(indicate)))
   if (a!="t")
     {plist$l <- c(plist$l,ind(A,1)[indicate] )
      plist$i <- c(plist$i,ind(A,2)[indicate] )
      plist$j <- c(plist$j,ind(A,3)[indicate])
     }
   else    # trend is a vector not an array
     {plist$l <- c(plist$l,rep(0,sum(indicate))) 
      plist$i <- c(plist$i, (1:length(A))[indicate]) 
      plist$j <- c(plist$j,rep(0,sum(indicate)))
     }
   plist  
 }# end locateARMA

       m <-dim(model$C)[3]
       p <-dim(model$A)[2]
       a <-dim(model$A)[1]
       b <-dim(model$B)[1]
       cc <-dim(model$C)[1]
       if (p!= dim(model$A)[3]) stop("Model A array dim 2 and 3 should be equal.")
       if (p!= dim(model$B)[2]) stop("Model B array dim inconsistent with array A.")
       if (p!= dim(model$B)[3]) stop("Model B array dim inconsistent with array A.")
       if (!is.null(model$C))
          if (p!= dim(model$C)[2]) stop("Model C array dim inconsistent with array A.")

       plist <- locateARMA(model$A,"A",p,p,a,
                     list(coef=NULL,location=NULL,i=NULL,j=NULL,
                         const=NULL,const.location=NULL,
                         const.i=NULL,const.j=NULL,l=NULL,const.l=NULL))
       plist <- locateARMA(model$B,"B",p,p,b,plist)
       if(!is.null(cc)) plist <- locateARMA(model$C,"C",p,m,cc,plist)
       if(!is.null(model$TREND)) 
            plist <- locateARMA(model$TREND,"t",p,m,cc,plist)

#       list.add(model, names(plist) ) <- plist
       model[ names(plist) ] <- plist
       model
} #end set.parameters.ARMA


set.arrays <- function(model, coefficients=NULL)  
 { # complete representaion info. based on parameter info.  
  UseMethod("set.arrays")
 }
   
set.arrays.TSestModel <- function(model, coefficients=NULL)  
 {set.arrays(TSmodel(model), coefficients=coefficients) }
    
set.arrays.SS <- function(model, coefficients=NULL){
	# N.B. Dimension and class (innov/ non.innov) info. is assumed accurate
    if (is.null(coefficients)) coefficients   <- coef(model)
    a.pos  <- model$location
    i.pos  <- model$i
    j.pos  <- model$j
    const  <- model$const
    ca.pos <- model$const.location
    ci.pos <- model$const.i
    cj.pos <- model$const.j  
    m <-dim(model$G)[2]
    n <-dim(model$F)[1]
    p <-dim(model$H)[1]
    f       <-  matrix(0,n,n)    # F     
    if(!is.null(m)) G        <-  matrix(0,n,m)    # control 
    H       <-  matrix(0,p,n)    # H       
    K        <-  matrix(0,n,p)    # Kalman gain
    if (is.null(model$Q)) Q <-  matrix(0,n,n)        # eventually discarded
    else                  Q <- array(0,dim(model$Q)) # system noise might not be nxn
    R        <-  diag(0,p,p)     # measurement noise
    z       <-  rep(0,n)        # initial state
    P       <-  diag(0,n)       # initial tracking error
    if (length(coefficients)>0) 
       {i <- a.pos == "f"
        f[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        if(!is.null(m)) 
          {i <- a.pos == "G"
           G[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
          }
        i <- a.pos == "H"
        H[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "K"
        K[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "Q"
        Q[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "R"
        R[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
        i <- a.pos == "z"
        z[i.pos[i]] <- coefficients[i]
        i <- a.pos == "P"
        P[cbind(i.pos[i],j.pos[i])] <- coefficients[i]
       }
    if (length(const)>0) 
       {i <- ca.pos == "f"
        f[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        if(!is.null(m)) 
          {i <- ca.pos == "G"
           G[cbind(ci.pos[i],cj.pos[i])] <- const[i]
          }
        i <- ca.pos == "H"
        H[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "K"
        K[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "Q"
        Q[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "R"
        R[cbind(ci.pos[i],cj.pos[i])] <- const[i]
        i <- ca.pos == "z"
        z[ci.pos[i]] <- const[i]
        i <- ca.pos == "P"
        P[cbind(ci.pos[i],cj.pos[i])] <- const[i]
       }
    model$F <- f 
    if(!is.null(m)) model$G <- G
    model$H <- H
    if (is.innov.SS(model))
      model$K <- K  
    else
      {model$Q <-Q
       model$R <-R
      }
    if(!is.null(model$z0)) model$z0<-z
    if(!is.null(model$P0)) model$P0<-P
    model 
} #end set.arrays.SS

set.arrays.ARMA <- function(model, coefficients=NULL) { 
	# N.B. Dimension and class info. is assumed accurate
       if (is.null(coefficients)) coefficients   <- coef(model)
       a.pos  <- model$location
       i.pos  <- model$i
       j.pos  <- model$j
       const  <- model$const
       ca.pos <- model$const.location
       ci.pos <- model$const.i
       cj.pos <- model$const.j  
       l.pos  <- model$l
       cl.pos <- model$const.l
       A  <-  array(0,dim(model$A))
       B  <-  array(0,dim(model$B))
       m <- dim(model$C)[3]
       p <- dim(model$A)[2]
       TREND <- rep(0,p)
       if(!is.null(m)) C  <-  array(0,dim(model$C))
       if (length(coefficients)>0) 
          {i <- a.pos == "A"
           A[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
           i <- a.pos == "B"
           B[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
           if(!is.null(m))
             {i <- a.pos == "C"
              C[cbind(l.pos[i],i.pos[i],j.pos[i])] <- coefficients[i]
             }
           i <- a.pos == "t"
           TREND[i.pos[i]] <- coefficients[i]
          }
    if (length(const)>0) 
          {i <- ca.pos == "A"
           A[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
           i <- ca.pos == "B"
           B[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
           if(!is.null(m))
             {i <- ca.pos == "C"
              C[cbind(cl.pos[i],ci.pos[i],cj.pos[i])] <- const[i]
             }
           i <- ca.pos == "t"
           TREND[ci.pos[i]] <- const[i]
          }
      model$A <- A
      model$B <- B
      if(!is.null(m)) model$C <- C 
      if(all(TREND==0)) model$TREND <- NULL
      else              model$TREND <-TREND
      model 
} #end set.arrays.ARMA


############################################################

# following is a workaround for ar in ts library

DSE.ar <- function(data, ...) {
  #fix for ar in R ts library (so that univariate case also gives array result)
  if (is.R()) if( !require("ts", warn.conflicts=FALSE)) stop("package ts is required.")
  res <- ar(ts(data), ...)
  if (! is.array(res$ar)) res$ar <- array(res$ar, c(length(res$ar),1,1))
  res
  }

print.test.value <- function(x, digits=16)
  {cat("c( ")
   if (all(is.na(x))) cat("NAs")
   else if (is.null(x)) cat("NULL")
   else if (is.logical(x)) cat(x)
   else if (is.S()) print(x, digits=digits)
   else 
     for (i in  1:length(x)) cat(", ", formatC(x[i], digits=digits, format="g"))
   cat(")\n")
   invisible()
  }

############################################################

#  Model simulation functions (to generate data)   <<<<<<<<<<

############################################################



simulate <- function(model, ...)UseMethod("simulate")

simulate.TSestModel <- function(model, input=input.data(model),
			sd=NULL, SIGMA=NULL, ...)
  {if (is.null(sd) & is.null(SIGMA)) SIGMA <- model$estimates$cov
   simulate(TSmodel(model), input=input, sd=sd, SIGMA=SIGMA, ...)
  }

simulate.SS <- function(model, input=NULL,
                 start. = NULL, freq = NULL,sampleT=100, 
                 noise=NULL, sd=1, SIGMA=NULL, rng=NULL, 
                 compiled=.DSECOMPILED)
{# S function to simulate a state space model:
#
#        z(t) = Fz(t-1) + Gu(t) + Qe(t)
#        y(t) = Hz(t)  + Rw(t)
# 
# or the innovations model:
#        z(t) = Fz(t-1) + Gu(t) + Kw(t-1)
#        y(t) = Hz(t)  + w(t)
#
#  see also the description in l.SS
# input=u must be specified if the matrix model$G is not NULL.
# If noise is NULL then an normal noise will be generated.
# This will be N(0,I) in the  non.innovation case (but Q and R 
# allow for arbitrary noise). If Q is not square (i.e. the system
# noise has a dimension less than the state dimension) then it is
# padded with zeros, so generated noise of higher dimension has no
# effect.  In the innovations case the noise will be N(0,sd^2).
# sd can be a vector of p elements corresponding to each of the p
# outputs.
# If noise is 
# specified it should be a list with elements $w0, $w and $e.
# $w0 is the noise at time zero (a p-vector of w(0) for innovations
# models and an n-vector of e(0) for  non.innovations models).
# If $w0 is a matrix (as for ARMA simulations) then it is set to a
# vector of zeros. This provides compatability with VAR models (ARMA
# models with no lags in B). In general ARMA and SS simulations will
# not produce exactly the same results because it is impossible to
# determine necessary transformation of initial conditions and w0.
# $w should be a sampleT by p matrix giving the output or 
# innovations noise or t=1 to sampleT. 
# For innovations models$e should be NULL.
# For  non.innovations models $e should be a sampleT by n matrix 
# giving the system noise for t=1 to sampleT.
# sampleT will be dim($w)[1] if noise is specified.

if (is.TSestModel(model)) model <- TSmodel(model)
if (!is.SS(model)) stop("simulate.SS expecting an SS TSmodel.")
 
 FF<-    model$F
 G <-    model$G
 H <-    model$H
 m <- dim(G)[2]
 if(is.null(m)) m <-0
 n <- dim(FF)[2]
 p <- dim(H)[1]
 if (m!=0)
   {if( is.null(input)) stop("input series must be supplied for this model.")
    if (sampleT != periods(input) ) input <- tfTruncate(input, end=sampleT)
   }
 if (is.innov.SS(model))
   {K <-    model$K}
 else
   {Q <-    model$Q
    if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square
    R <-    model$R
   }
 
 if(is.null(rng)) rng <- set.RNG() # returns setting so don't skip if NULL
 else        {old.rng <- set.RNG(rng);  on.exit(set.RNG(old.rng))  }
 

set.ts <- TRUE             
if (!is.null(start.))
  {if (!is.null(freq))   tf <- list(start=start., frequency=freq)
   else
      {warning("start. is specified but not freq. Using freq=1.") 
       tf <- list(start=start., frequency=1)
      }
  }  
else if( (!is.null(input))   && is.tframed(input))   tf <- tframe(input)
else if ((!is.null(noise$w)) && is.tframed(noise$w)) tf <- tframe(noise$w)
else set.ts <-  FALSE



# It would be better to use makeTSnoise here (lag=1 and w0<-c(noise$w0) ) and
#  possibly two calls to get e for non-innov models. However, this will affect
#  historical comparisons so it needs to be done carefully!
 e <-NULL
 if (is.null(noise))
   {if (!is.innov.SS(model)) 
      {w0 <- rnorm(n)
       e <- matrix(rnorm(sampleT*n),sampleT,n)
       w <- matrix(rnorm(sampleT*p),sampleT,p)
      }
    else 
      {w <- matrix(NA,sampleT,p)
       w0 <- rep(NA,p)
       if (!is.null(SIGMA))
         {if(length(SIGMA) == 1) SIGMA <- diag(SIGMA, p)
	  W <- t(chol(SIGMA))
	  w.help <- t(W %*% matrix(rnorm((sampleT+1)*p),nrow=p,ncol=sampleT+1))
	  w0 <- w.help[1,]
	  w <- w.help[-1,]
         }
       else
         {if (length(sd)==1) sd <-rep(sd,p)
          for (i in 1:p)
            {w0[i] <- rnorm(1,sd=sd[i])
             w[,i] <- rnorm(sampleT,sd=sd[i])
            }
         }
      }
   }
 else
   {w0 <- noise$w0
    if (is.matrix(w0)) w0 <- rep(0,p) # see note above re VAR
    w<-noise$w
    e<-noise$e
    sampleT<-dim(noise$w)[1]
   }

 y <- matrix(0,sampleT,p)
 state <- matrix(0,sampleT,n)
 if(is.null(model$z0)) z <-rep(0,n) # initial state
 else  z <- model$z0        

 if (is.innov.SS(model)) 
   {z <-  c(FF%*% z) + c(K %*% w0)
    if (m !=0) z <-  z + c(G%*%input[1,]) 
    y[1,]  <-  c(H %*% z) + c(w[1,])
   }   
 else      
   {z <-  c(FF%*% z) + c(Q %*% w0)
    if (m !=0) z <-  z + c(G %*% input[1,])  
    y[1,]  <-  c(H %*% z) + c(R %*% w[1,])
   }
 state[1,]<-z
 # Note: the first period is done above for both compiled and S versions so
 #    that initial conditions do not have to be passed to compiled code.
 if (compiled)
   {if (is.innov.SS(model))
      {Q <-matrix(0,n,n)
       R <-matrix(0,p,p)
       e <- matrix(0,sampleT,n)
      }
    else K <- matrix(0,n,p)
    if (m==0) 
      {input <- matrix(0,sampleT,1)
       G <- matrix(0,n,1)
      }
    storage.mode(y)     <- "double"
    storage.mode(state) <- "double"
    r<- .Fortran("simss",y=y, 
                         state=state, 
                         as.integer(m),
                         as.integer(n),
                         as.integer(p), 
                         as.integer(sampleT),  
                         as.double(input),
                         as.double(w),
                         as.double(e),
                         as.double(FF),
                         as.double(G),   
                         as.double(H),
                         as.double(K), 
                         as.double(Q),      
                         as.double(R),    
                         as.integer(is.innov.SS(model)), DUP=.DSEDUP
			 )[c("y","state")]
    y <- r$y
    state <- r$state
    if (m==0) input <- NULL
   }
 else
   {for (Time in 2:sampleT)  
      {z <-  c(FF%*% z) 
       if (m !=0) z <-  z + c(G%*%input[Time,]) 
       if (is.innov.SS(model))  z <-  z + c(K%*%w[Time-1,])
       else       z <-  z + c(Q%*%e[Time-1,])
       state[Time,] <- z
       if (is.innov.SS(model)) y[Time,]  <-  c(H %*% z) + c(w[Time,])   
       else      y[Time,]  <-  c(H %*% z) + c(R %*% w[Time,])
   }  }
 if (set.ts)
   { y     <- tframed(y,     tf=tf, names=seriesNamesOutput(model)) 
     state <- tframed(state, tf=tf )
   }
 else  seriesNames(y) <- seriesNamesOutput(model)
 TSdata(list(input=input,output=y, state=state, version=version, 
   model=model, description="data generated by simulate.ss",
   noise=list(w0=w0,w=w, e=e, rng=rng, SIGMA=SIGMA, sd=sd)))
}

simulate.ARMA <- function(model, y0=NULL, input=NULL, input0=NULL,
                start. = NULL, freq = NULL, sampleT=100,
                noise=NULL, sd=1, SIGMA=NULL,
                rng=NULL, noise.model=NULL, 
                compiled=.DSECOMPILED)
{# S function to simulate ARMA mode:
#
#       A(L)y(t) =  B(L)w(t) + C(L)u(t) + TREND
# 
# See also the description in ARMA.s
# input=u must be specified if the matrix model$C is not NULL.
# The rng will be set first if it is specified. If noise is 
# specified it should be a list with elements $w0 and $w.
# $w0 is the noise(w) prior to time=1 (a (dim($B)[1]-1) by p matrix).
# For VAR models B has no lags so w0 has no effect.
# $w should be a sampleT by p matrix giving the 
# noise for t=1 to sampleT. 
# sampleT will be dim($w)[1] if noise is specified.

if (is.TSestModel(model)) model <- TSmodel(model)
if (!is.ARMA(model)) stop("simulate.ARMA expecting an ARMA TSmodel.")
 
A<-    model$A
B <-    model$B
C <-    model$C
TREND <- model$TREND
m <- dim(C)[3]
if (is.null(m)) m <-0
p <- dim(A)[2]
a <-dim(A)[1]
b <-dim(B)[1]
if (is.null(C)) cc <- 0
else            cc<-dim(C)[1]
if ( (p != dim(A)[3])
    |(p != dim(B)[2])
    |(p != dim(B)[3])) 
      stop("dimension of model parameters do not conform!")
if (0 !=m) if (p != dim(C)[2]) 
      stop("dimension of model parameters do not conform!")
if (m!=0)
   {if( is.null(input)) stop("input series must be supplied for this model.")
    if (sampleT != periods(input) ) input <- tfTruncate(input, end=sampleT)
   }
 
if (p==1) invA0 <- matrix(1/A[1,,],1,1)
else      invA0 <- solve(A[1,,])
for (l in 1:a) A[l,,] <- invA0 %*% A[l,,]      # set A(0) = I      
for (l in 1:b) B[l,,] <- invA0 %*% B[l,,] 
if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invA0 %*% C[l,,]  
if(!is.null(TREND)) TREND <- invA0 %*% TREND

set.ts <- TRUE             
if (!is.null(start.))
  {if (!is.null(freq))   tf <- list(start=start., frequency=freq)
   else
      {warning("start. is specified but not freq. Using freq=1.") 
       tf <- list(start=start., frequency=1)
      }
  }  
else if( (!is.null(input))   && is.tframed(input))   tf <- tframe(input)
else if ((!is.null(noise$w)) && is.tframed(noise$w)) tf <- tframe(noise$w)
else set.ts <-  FALSE

noise <- makeTSnoise(sampleT,p,b, noise=noise, rng=rng,
                        SIGMA=SIGMA, sd=sd, noise.model=noise.model,
                        start=start., frequency=freq)
if (is.null(sampleT)) sampleT<-noise$sampleT
 
 if(is.null(y0)) y0<-matrix(0,a,p)
 if((m!=0) & is.null(input0)) input0 <- matrix(0,dim(C)[1],m)
 y <- matrix(0,sampleT,p)  
 if (compiled)
   {if (m==0) 
      {input <- matrix(0,sampleT,1)
       input0 <- matrix(0,1,1)
       C <- matrix(0,1,1)
      }
    if (is.null(TREND)) TREND<- rep(0,p)
#    yo<- list(y=y, y0,m,p, a, b, cc, sampleT,input,input0,w,w0,A,B, C,TREND)
    storage.mode(y)     <- "double"

    y <-.Fortran("simrma",y=y, 
                         as.double(y0),
                         as.integer(m),
                         as.integer(p), 
                         as.integer(a), 
                         as.integer(b), 
                         as.integer(cc), 
                         as.integer(sampleT),  
                         as.double(input[1:sampleT,]),
                         as.double(input0),
                         as.double(noise$w),
                         as.double(noise$w0),
                         as.double(A),
                         as.double(B),   
                         as.double(C),
                         as.double(TREND), DUP=.DSEDUP) [["y"]]

    if (m==0) 
      {input  <- NULL
       input0 <- NULL
      }
   }
 else
  {w0 <- noise$w0
   w<-noise$w
   for (Time in 1:sampleT)  
   {if(!is.null(TREND)) y[Time,] <- TREND # + y[Time,] 
    for (l in 2:a) 
       if(Time+1-l<=0)
          if (p==1) y[Time,] <- y[Time,]-c(A[l,,]  *  y0[l-Time,]) 
          else      y[Time,] <- y[Time,]-c(A[l,,] %*% y0[l-Time,])
       else                
          if (p==1) y[Time,] <- y[Time,]-c(A[l,,]  *  y[Time+1-l,])
          else      y[Time,] <- y[Time,]-c(A[l,,] %*% y[Time+1-l,])
    for (l in 1:b) 
       if (Time+1-l<=0)
          if (p==1) y[Time,]<- y[Time,] +c(B[l,,]  *  w0[l-Time,]) 
          else      y[Time,]<- y[Time,] +c(B[l,,] %*% w0[l-Time,])
       else
          if (p==1) y[Time,]<- y[Time,] +c(B[l,,]  *  w[Time+1-l,])
          else      y[Time,]<- y[Time,] +c(B[l,,] %*% w[Time+1-l,])
    if (m!=0) for (l in 1:cc) 
       if (Time+1-l<=0)
          if (m==1) y[Time,]<- y[Time,] + c(C[l,,]  *  input0[l-Time,])
          else      y[Time,]<- y[Time,] + c(C[l,,] %*% input0[l-Time,])
       else
          if (m==1) y[Time,]<-y[Time,] + c(C[l,,]  *  input[Time+1-l,])
          else      y[Time,]<-y[Time,] + c(C[l,,] %*% input[Time+1-l,])
  }}
 if (set.ts) y <- tframed(y, tf=tf, names=seriesNamesOutput(model)) 
 else seriesNames(y) <- seriesNamesOutput(model)
 TSdata(list(input=input,output=y, 
          model=model, input0=input0, 
          description="data generated by simulate.ARMA", 
#          prior.args=prior.args, post.args=post.args,
          noise=noise))
}

#######################################################################

#functions which work on data (and models)  <<<<<<<<<<

############################################################

#     likelihood and residual calculation functions  <<<<<<<<<<

############################################################


#  L <- function(residual)
#  { # negative log likelihood of a residual
#    sampleT <-nrow(residual)
#    p <- ncol(residual)
#  #  Om <- var(residual)  # var removes mean and /sampleT-1
#    Om <-t(residual) %*% residual /sampleT
#    v <- svd(Om) #eigenvalues are not robust to degenerate density.
#    # next seem a hard way to get det
#  # like1 <- 0.5 * sampleT * log(prod(v$d[v$d!=0]))
#    like1 <- 0.5 * sampleT * log(prod(
#                 v$d[v$d > (v$d[1]*sqrt(.Machine$double.eps))]))
#  # svd is more robust than solve(Om) for degenerate densities
#  #  if (1 == length(v$d)) OmInv <-  v$v %*% (1/v$d) %*%t(v$u) 
#  #  else OmInv <-  v$v %*% diag(1/v$d) %*%t(v$u) 
#  #  OmInv <-  v$v %*% (t(v$u) * 1/v$d) is faster and equivalent but relies on 
#  #  recycling of d and columnwise storage, which work in S and R but are "tricks"
#  # A better way is (but should account for degenerate space as in like1 above)
#    OmInv <-  v$v %*% sweep(t(v$u),1,1/v$d, "*") 
#  #  like2 <- sum(diag(residual %*%OmInv %*% t(residual))) /2
#    like2 <- sum(residual * (residual %*% OmInv)) /2
#    const <- (sampleT * p * log(2 * pi))/2
#    c(const+like1+like2, const, like1,like2)
#  }

residual.stats <- function(pred, data, sampleT=nrow(pred), warn=TRUE)
{  # pred and data should be matrices (model prediction and output data).
   # sampleT allows for the possibility that a sub-sample of data 
   #   was used for estimation.
   # (note predictT can be determined from nrow(pred) and is not used.)
   e <- if (is.null(pred))     -data[1:sampleT,,drop=FALSE]
        else if (is.null(data)) pred[1:sampleT,,drop=FALSE]
        else               pred[1:sampleT,,drop=FALSE] - data[1:sampleT,,drop=FALSE]
   p <- ncol(e)
   Om <-t(e) %*% e /sampleT
   if (any(is.na(Om))) {like1 <- like2 <- 1e100}  
   else if (any(Om >1e100)) {like1 <- like2 <- 1e100}  
   else
     {v <- svd(Om) #eigenvalues are not robust to degenerate density.
#     i <- v$d!=0
     i <- v$d > (v$d[1]*sqrt(.Machine$double.eps))
      if (!all(i))
        {if(warn) warning("The cov. matrix is singular. Working on subspace.")
         v$d <- v$d[i]
         v$u <- v$u[i,i, drop=FALSE]
         v$v <- v$v[i,i, drop=FALSE]
         e <- e[,i, drop=FALSE]
        }
      like1 <- 0.5 * sampleT * log(prod(v$d)) # det
#      if (1 == length(v$d)) OmInv <-  v$v %*% (1/v$d) %*%t(v$u) 
#      else OmInv <-  v$v %*% diag(1/v$d) %*%t(v$u) # more robust than solve(Om)
#  following is equivalent
#      OmInv <-  v$v %*% sweep(t(v$u),1,1/v$d, "*") 
#      like2 <- sum(e * (e %*% OmInv)) /2
#  but this works out to (sampleT*p/2) and fixing for degenerate distributions:
       like2 <- (sampleT*length(v$d))/2
     }
   const <- (sampleT * p * log(2 * pi))/2
   invisible(list(like=c(const+like1+like2, const, like1,like2),
                  cov =Om, pred=pred, sampleT=sampleT))
}



sum.sqerror <- function(coefficients, model=NULL, data=NULL, error.weights=NULL) 
{ #  this returns only the sum of the weighted squared errors (eg.for optimization).
#  If model, data or error.weights are not supplied the program looks for
#    a global variable named Obj.Func.ARGS with corresponding elements.
#  The sample size is determined by periodsOutput(data).
 if ( is.null(model)) stop("model missing") # model <- Obj.Func.ARGS$model
 if ( is.null(data))  stop("data missing") # data  <- Obj.Func.ARGS$data
 if ( is.null(error.weights)) stop("error.weights missing") #error.weights <- Obj.Func.ARGS$error.weights 
 sum(l(set.arrays(model,coefficients=coefficients), data,
       result="weighted.sqerror",error.weights=error.weights))
}




l <- function(obj1, obj2, ...)UseMethod("l")
l.TSdata <- function(obj1, obj2, ...) {l(obj2, obj1, ...) }
l.TSestModel <- function(obj1, obj2, ...) {l(TSmodel(obj1), obj2, ...)}


l.ARMA <- function(obj1, dat, sampleT=NULL, predictT=NULL,result=NULL,
                error.weights=0,  compiled=.DSECOMPILED, 
		warn=TRUE, return.debug.info=FALSE)
{#  calculate likelihood, residuals, prediction, etc. for ARMA model
 # N.B.  The compiled version is much preferred for speed.
 #  sampleT is the length of data which should be used for 
 #  calculate the one-step ahead predictions, and likelihood value for the model:
#
#       A(L)y(t) =  B(L)w(t) + C(L)u(t) + TREND
# 
# A(L) (axpxp) is the auto-regressive polynomial array.
# B(L) (bxpxp) is the moving-average polynomial array.
# C(L) (cxpxm) is the  input polynomial array.
# TREND is a constant vector added at each period.
# y is the p dimensional output data.
# u is the m dimensional control (input) data.
# Om is the estimated output cov matrix.

model <- if (is.TSestModel(obj1)) TSmodel(obj1) else obj1
if (!is.ARMA(model)) stop("l.ARMA expecting an ARMA TSmodel.")

dat <- freeze(dat)
if(!check.consistent.dimensions(model,dat)) stop("dimension error")
if (is.null(sampleT))  sampleT  <- periodsOutput(dat)
if (is.null(predictT)) predictT <- sampleT
if (sampleT > predictT) stop("sampleT cannot exceed predictT.")
if (sampleT > periodsOutput(dat)) stop("sampleT cannot exceed length of data.")
if (0 != nseriesInput(dat))
  {if (periodsInput(dat) < predictT)
      stop("input data must be at least as long as requested prediction.")
   if (any(is.na(input.data(dat)))) stop("input data cannot contain NAs.")
  }
if (any(is.na(output.data(dat)))) stop("output data cannot contain NAs.")

u <- input.data(dat)
y <- output.data(dat)
A<-    model$A
B <-   model$B
C <-   model$C
TREND <- model$TREND
m <- dim(C)[3]
if (is.null(m)) m <-0
p <- dim(A)[2]
a <-dim(A)[1]
b <-dim(B)[1]
if ( (p != dim(A)[3])
    |(p != dim(B)[2])
    |(p != dim(B)[3])) 
      stop("dimension of model parameters do not conform!")
if (0 !=m) if (p != dim(C)[2]) 
      stop("dimension of model parameters do not conform!")
if (p != nseriesOutput(dat))
      stop("dimension of model parameters do not conform with the data!")
if (m == 0) 
   {if(!is.null(u)) 
      stop("model parameters indicate an no input but input data exists!")
   }
else if (m != dim(u)[2])  
      stop("dimension of model parameters do not conform with the data!")
if (compiled)
  {if (m==0)
     {C <- array(0,c(1,p,1))    # can't pass 0 length array to compiled
      u <- matrix(0,predictT,1)
     }
   if (is.null(model$TREND)) TREND <- rep(0,p)
   is  <- max(m,p)
   r  <- .Fortran("arma",
                         pred=matrix(1e20,predictT,p), # pred,     
                         as.integer(length(error.weights)),
                         weighted.sqerror=matrix(0,sampleT,p),
                         error.weights= as.double(error.weights),
                         as.integer( m), 
                         as.integer( p) ,      
                         as.integer( dim(A)[1]),  # 1+order of A  
                         as.integer( dim(B)[1]),  # 1+order of B  
                         as.integer( dim(C)[1]),  # 1+order of C  
                         sampleT=as.integer(sampleT),
                         predictT=as.integer(predictT),
                         as.integer(periodsOutput(dat)), 
                         as.double(u), # as.double() works ok with compiled but
                                    #messes up the dim(u) returned in the list
                         as.double(y),         
                         as.double(A),  
                         as.double(B),   
                         as.double(C),
                         as.double(TREND),
                         as.integer(is),  # scratch array dim
                         as.double(matrix(0,is,is)),  # scratch array
                         as.double(matrix(0,is,is)),  # scratch array
                         as.double(rep(0,is)),         # scratch array
                         DUP=.DSEDUP) [c("pred", "weighted.sqerror")]
   if (all(0==error.weights)) r$weighted.sqerror <- NULL
  }
else   # start S version
  {prederror <- matrix(0,sampleT,p)  
   wt.err <- NULL
   invB0 <- solve(B[1,,])
   for (l in 1:a) A[l,,] <- invB0 %*% A[l,,]      # set B(0) = I      
   for (l in 1:b) B[l,,] <- invB0 %*% B[l,,]  
   if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invB0 %*% C[l,,]  
   if(!is.null(TREND)) TREND <- invB0 %*% TREND
   if (1 < length(error.weights)) wt.err <- matrix(0,predictT,p)    
   for (Time in 1:sampleT)  
      {if(!is.null(TREND)) vt <- -TREND
       else vt    <-  rep(0,p) 
       for (l in 1:a)
          if (l<=Time)  #this is cumbersome but drop=FALSE leaves A,B,C as 3 dim arrays
              if (p==1) vt <- vt + c(A[l,,]  *  y[Time+1-l,])
              else      vt <- vt + c(A[l,,] %*% y[Time+1-l,])  
       if (b >= 2) for (l in 2:b) 
          if (l<=Time) 
             if (p==1) vt <- vt - c(B[l,,]  *  prederror[Time+1-l,]) 
             else      vt <- vt - c(B[l,,] %*% prederror[Time+1-l,])
       if (m!=0) for (l in 1:dim(C)[1]) 
          if (l<=Time) 
             if (m==1) vt <- vt - c(C[l,,]  *  u[Time+1-l,])
             else      vt <- vt - c(C[l,,] %*% u[Time+1-l,])  
       prederror[Time,] <- vt #this is not really the pred error unless B0=I   

       if (any(0!=error.weights))          
        {wt.err[Time,] <- error.weights[1]*vt^2  # weighted sq prediction error
         if (length(error.weights)>1)
           {for (h in 2:length(error.weights))
            if ( (Time+h-1) <= sampleT)
              {if(!is.null(TREND)) vt <- -TREND
               else vt    <-  rep(0,p) 
               for (l in 1:a)
                  if (l < Time+h) 
                     if (p==1) vt <- vt + c(A[l,,]  *  y[Time+h-l,])
                     else      vt <- vt + c(A[l,,] %*% y[Time+h-l,])  
               if (b >= 2) for (l in 2:b) 
                  if (l < Time+h) 
                     if (p==1) vt <- vt - c(B[l,,]  *  prederror[Time+h-l,]) 
                     else      vt <- vt - c(B[l,,] %*% prederror[Time+h-l,])
               if (m!=0) for (l in 1:dim(C)[1]) 
                   if (l < Time+h) 
                      if (m==1) vt <- vt - c(C[l,,]  *  u[Time+h-l,])
                      else      vt <- vt - c(C[l,,] %*% u[Time+h-l,]) 
               wt.err[Time,] <- wt.err[Time,] + error.weights[h]*(solve(invB0 )%*%vt)^2
     }  }  }
  }

   pred <- matrix(0,predictT,p)
   pred[1:sampleT,] <- y[1:sampleT,,drop=FALSE] - prederror[1:sampleT,,drop=FALSE] %*% t(solve(invB0)) 
   # now multi-step predictions to predictT
   if (predictT > sampleT)
    {for (Time in (sampleT+1):predictT)  
      {invA0 <- solve(A[1,,])
       for (l in 1:a) A[l,,] <- invA0 %*% A[l,,]      # set A(0) = I      
       for (l in 1:b) B[l,,] <- invA0 %*% B[l,,]  
       if (m!=0) for (l in 1:dim(C)[1]) C[l,,] <- invA0 %*% C[l,,]  
       if(!is.null(TREND)) TREND <- invA0 %*% TREND
       if(!is.null(TREND)) pred[Time,] <- pred[Time,]+TREND
       for (l in 2:a) 
          if(Time+1-l<=sampleT)
             if (p==1) pred[Time,] <- pred[Time,]-c(A[l,,]  *  y[Time+1-l,]) 
             else      pred[Time,] <- pred[Time,]-c(A[l,,] %*% y[Time+1-l,])
          else                
             if (p==1) pred[Time,] <- pred[Time,]-c(A[l,,]  *  pred[Time+1-l,])
             else      pred[Time,] <- pred[Time,]-c(A[l,,] %*% pred[Time+1-l,])
       if (b >= 2) for (l in 2:b) 
          if (Time+1-l<=sampleT)
             if (p==1)  pred[Time,] <- pred[Time,] +c(B[l,,]  *  prederror[Time+1-l,]) 
             else       pred[Time,] <- pred[Time,] +c(B[l,,] %*% prederror[Time+1-l,])
       if (m!=0) for (l in 1:dim(C)[1]) 
          if (m==1) pred[Time,] <- pred[Time,] + c(C[l,,]  *  u[Time+1-l,]) 
          else      pred[Time,] <- pred[Time,] + c(C[l,,] %*% u[Time+1-l,])
      }
    }
   r<-list(pred=pred,   weighted.sqerror=wt.err)
  } # end of S version

tf <- tfTruncate(tframe(output.data(dat)), end=predictT)
tframe(r$pred) <- tf
if (! is.null(r$weighted.sqerror))  tframe(r$weighted.sqerror) <- tf
if((!is.null(result)) && (result == "pred")) return(r$pred)
r <- append(residual.stats(r$pred, y, sampleT, warn=warn), 
        list(error.weights=error.weights, weighted.sqerror=r$weighted.sqerror))

if (return.debug.info) 
    r$debug.info <-list(m=m, p=p, a=a,b=b,c=c, A=A,B=B,C=C,TREND=TREND, 
          u=u, y=y,
          prederror=prederror, pred=pred, invB0=invB0, wt.err=wt.err, 
          error.weights=error.weights, sampleT=sampleT)

if ( is.null(result)) return(classed( # TSestModel constructor (l.ARMA)
              list(estimates=r, data=dat, model=model), "TSestModel"))
else 
   {if (result =="like") return(r$like[1]) # neg.log.like. from residual.stats
    else { return(r[[result]]) }
   }
stop("should never get to here in l.ARMA.")
}                       


l.SS <- function(obj1, data, sampleT=NULL, predictT=NULL, error.weights=0,
                 return.state=FALSE, return.track=FALSE, result=NULL, 
		 compiled=.DSECOMPILED,
                 warn=TRUE, return.debug.info=FALSE)
{# ref. B.D.O.Anderson & J.B.Moore "Optimal Filtering" p.39,44.
# sampleT is the length of data which should be used for calculating
# one step ahead predictions. y must be at least as
#  long as sampleT. If predictT is large than sampleT then the model is simulated to 
# predictT. y is used if it is long enough. u must be at least as long as predictT.
# The default result=0 returns a list of all the results. Otherwise only the 
#    indicated list element is return (eg. result=1 return the likelihood and
#    result=3 returns the one step ahead predictions.
#   Calculate the state, residuals, and likelihood value for the model:
#
#        z(t) = Fz(t-1) + Gu(t) + Qe(t)
#        y(t) = Hz(t)  + Rw(t)
# 
# or the innovations model:
#        z(t) = Fz(t-1) + Gu(t) + Kw(t-1)
#        y(t) = Hz(t)  + w(t)
#
# FF (nxn) is the state transition matrix F.
# H (pxn)is the output matrix H.
# Q (nxn) is the input matrix of the system noise and the noise is assumed to be white. 
#    Some authors (eg. Harvey) modify this as rt*qt*rt' where rt is the matrix for the 
#    system noise and qt is the noise cov, but that is redundant.
# R (pxp) is the input matrix of the output (measurement) noise, which is assumed white. 
#      probably need R if p>n ??
# G (nxp)is the control (input) matrix.
# K (nxp)is the Kalman gain.
# y is the p dimensional output data.
# u is the m dimensional exogenous (input) data.
# z is the n dimensional (estimated) state at time t,  E[z(t)|y(t-1), u(t)] denoted E[z(t)|t-1].
#    Note: In the case where there is no input u this corresponds to what
#     would usually be called the predicted state - not the filtered state.
# state is the history of the state.
# Om is the estimated output cov matrix.
# vt is the prediction error.
# pred is the history of the one-step ahead predictions, E[y(t)|y(t-1),u(t)] denoted E[y(t)|t-1].
# The history of the prediction error is given by y-pred[1:predictT,]or y-pred[1:sampleT,]
#     If error.weights is greater than zero then weighted prediction 
#     errors are calculated up to the horizon indicated
#     by the length of error.weights. The weights are applied to the squared
#     error at each period ahead.
# P is the one step ahead estimate of the state tracking error matrix at each 
# period. Cov{z(t)-E[z(t)|t-1]}
# trackError is the history of P.
#       Tracking error pt can only be calculated if Q and R are provided ( gain FALSE).
#       Using the Kalman gain K directly these are not necessary 
#       for the likelihood calculation,
#       but the tracking error cannot be calculated.
# If z0 is supplied it is used as the estimate of the state at time 0.
# If not supplied it is set to zero.
# If P0 is supplied it is used as the initial tracking error P(t=1|t=0).
# If not supplied it is set to I.
# could check that Q is symmetric  or positive definite but ...

model <- if (is.TSestModel(obj1)) TSmodel(obj1) else obj1
if (!is.SS(model)) stop("l.SS expecting an SS TSmodel.")

data <- freeze(data)
if(!check.consistent.dimensions(model,data)) stop("dimension error\n")
if (is.null(sampleT))  sampleT  <- periodsOutput(data)
if (is.null(predictT)) predictT <- sampleT
if (sampleT > predictT) stop("sampleT cannot exceed predictT.\n")
if (sampleT > periodsOutput(data))
    stop("sampleT cannot exceed length of data.\n")
if (0 != nseriesInput(data))
  {if (periodsInput(data) < predictT)
      stop("input data must be at least as long as requested prediction.\n")
   if (any(is.na(input.data(data)))) stop("input data cannot contain NAs.\n")
  }
if (any(is.na(output.data(data)))) stop("output data cannot contain NAs.\n")

gain <- is.innov.SS(model)
if (gain & return.track) 
   warning("Tracking error is zero for an innovations model. track will not be calculated.")

FF<-    model$F
H <-    model$H
n <- dim(FF)[2]
p <- dim(H)[1]
if (is.null(model$G))
  {m<-0
   G<-matrix(0,n,1)       # can't call compiled with 0 length arrays
   u <- matrix(0,predictT,1)
  }
else
  {m <- dim(model$G)[2]
   G <-model$G
   u <- input.data(data)
  } 
if (gain)            # K or Q,R can be NUll in model, which messes up compiled
   {K <-    model$K
    Q <-    matrix(0,1,1)      #not used
    R <-    matrix(0,1,1)      #not used
    track <-array(0,c(1,1,1))  #not used
   }
else
   {Q <-    model$Q
    if (ncol(Q)<n) Q <- cbind(Q,matrix(0,n,n-ncol(Q))) # Q assumed square in compiled
    R <-    model$R
    K <-    matrix(0,n,p)      # this is used
    if(return.track) track <-array(0,c(predictT,n,n))
    else             track <-array(0,c(1,1,1))  #not used
    storage.mode(track) <- "double"
   }
if (return.state | return.debug.info) state <- matrix(0,predictT,n)
else              state <- matrix(0,1,1)    #not used
storage.mode(state) <- "double"
if(is.null(model$z0)) z <-rep(0,n)   # initial state
else  z <-model$z0
if(is.null(model$P0)) P <- diag(1,n) # initial state tracking error 
else  P <-model$P0               # this is not used in innov. models

if (compiled)
  {r <- .Fortran("kf",
                  pred=matrix(0,predictT,p),    
                  as.integer(length(error.weights)), 
                  weighted.sqerror=matrix(0,sampleT,p),
                  error.weights=as.double(error.weights),   
                  as.integer(return.state),
                  state=state,         
                  as.integer(return.track & !gain),
                  track=track,                  
                  as.integer(m), 
                  as.integer(n), 
                  as.integer(p), 
                  sampleT=as.integer(sampleT), 
                  predictT=as.integer(predictT), 
                  as.integer(periodsOutput(data)),  
                  as.double(u), 
                  as.double(output.data(data)),  
                  as.double(FF),   
                  as.double(G),   
                  as.double(H),  
                  as.double(K), 
                  as.double(Q),      
                  as.double(R),    
                  as.integer(gain),
                  as.double(z),
                  as.double(P), 
		  DUP=.DSEDUP) [c("pred","state","track","weighted.sqerror")]
   if (all(0==error.weights)) r$weighted.sqerror <- NULL
  }
else                  #  S version
  {y <- output.data(data)
   vt    <-  rep(0,p)     # initial prediction error
   pred  <- matrix(0,predictT,p) 
   wt.err <- NULL
   if (1 < length(error.weights)) wt.err <- matrix(0,predictT,p)
   if ( ! gain ) 
     {RR <- R %*% t(R)  
      QQ <- Q %*% t(Q)  
     }                            
                                   
   for (Time in 1:sampleT)  {
       if ( ! gain) 
         {PH  <-  P %*% t(H)
          ft    <- ( H %*% PH )  + RR         
          ft    <-  (ft + t(ft))/2   # force ft to be symmetric 
          K   <-  t(solve(ft,t(FF %*% PH)))  
          P   <-  (FF %*% P %*% t(FF) ) - ( K %*% H %*% P %*% t(FF) ) + QQ  # P(t|t-1)
          P   <-  (P + t(P))/2  # force symmetry (eliminate rounding error problems)
          if (return.track) track[Time,,] <- P   # P(t|t-1)
          # note P(t|t) = P-P%*%t(H)%*%solve(H%*/home/mfa5/gilp/dse/my/src/SCCS/s.dse1b.hs*%t(H)+RR)%*04/28/00*%P
         } 
         
       z<- c(FF%*%z) + c(K%*%vt)  # E[z(t)| t-1 ]
       if (m !=0) z<- z + c(G%*%u[Time,])
       if (return.state | return.debug.info) state[Time,]<- z
       pred[Time,] <- Ey  <-  c(H %*% z)       # predicted output     
       vt<-  y[Time,] - Ey                     # prediction error 
       if (any(0!=error.weights))          
        {wt.err[Time,] <- error.weights[1]*vt^2  # weighted sq prediction error
         if (length(error.weights)>1)
          {zh <-z
           for (h in 2:length(error.weights))
            if ( (Time+h-1) <= sampleT)
              {zh <-  c(FF%*%zh)
               if (h==2) zh <- zh + c(K%*%vt) # vt is 0 for h>2
               if (m !=0) zh<- zh + c(G%*%u[Time+h-1,])
               wt.err[Time,] <- wt.err[Time,] + 
                          error.weights[h]*(y[Time+h-1,] -  c(H %*% zh))^2
   }   }  }   }


#   prederror <- y[1:sampleT,,drop=FALSE]-pred[1:sampleT,,drop=FALSE]

   # now multi-step prediction to predictT 
   # This requires u but not y (y is ignored if it is supplied)
   if (predictT > sampleT)
    {for (Time in (sampleT+1):predictT)  
       {z <-  c(FF%*% z) 
        if (m !=0) z <-  z + c(G%*%u[Time,]) 
        if (Time==sampleT+1) z <- z + c(K%*%vt)
        if (return.state) state[Time,] <- z
        pred[Time,]  <-  c(H %*% z)                  # predicted output 
       }
     }
   r<- list(pred=pred, state=state, track=track, weighted.sqerror=wt.err) 
  }   # end of S version

tf <- tfTruncate(tframe(output.data(data)), end=predictT)
tframe(r$pred) <- tf
if (! is.null(r$weighted.sqerror))  tframe(r$weighted.sqerror) <- tf

filter <-NULL
if (return.state | return.track)
  {if (gain|(!return.track))  filter$track <- NULL 
   else                       
     {filter$track <- r$track
      tframe(filter$track) <- tf
     }
   if (return.state)
     {filter$state <- r$state
      tframe(filter$state) <- tf
     }
  }
if((!is.null(result)) && (result == "pred")) return(r$pred)
r <- append(residual.stats(r$pred, output.data(data), sampleT, warn=warn), 
        list(error.weights=error.weights, weighted.sqerror=r$weighted.sqerror))

if (return.debug.info) 
   r$debug.info <- list(m=m, p=p, a=a,b=b,c=c, G=G,FF=FF,K=K,P=P, H=H, u=u,y=y,
        pred=pred, wt.err=wt.err, error.weights=error.weights, sampleT=sampleT)

if ( is.null(result)) 
   {r <-list(estimates=r, data=data, model=model, filter=filter) 
    return( classed(r, "TSestModel")) # constructor (l.SS)
   }
else 
   {if (result =="like") return(r$like[1]) # neg.log.like. from residual.stats
    else { return(r[[result]]) }
   }
"should never get to here in l.SS"
} # end of l.SS


smoother <- function(model, data, compiled=.DSECOMPILED){
 #  Fixed interval smoother for a model as returned by l.SS.
 # ref. appendix of Shumway and Stoffer,1982, J.of Time Series, 253-264,
 #        Jazwinski 1970, or Anderson and Moore.
 # Note: this does not allow the same option as l.SS for calculating over a
 #    sub-sample. Smoothing is done over the length of the available filter
 #    data (which will be calculated to the length of the data if not
 #    supplied). For models with an input smoothing will only be done to the
 #    length of input data if that is smaller than the available filter data. 
 # See l.SS for details of the model:
 #
 #        z(t) = Fz(t-1) + Gu(t) + Qe(t)
 #        y(t) = Hz(t)  + Rw(t)
 # 
 filter <- NULL
 estimates <- NULL
 if (is.TSestModel(model)) 
   {filter    <- model$filter
    estimates <- model$estimates
    model     <- TSmodel(model)
   }
 if  (!is.non.innov.SS(model)) stop("smoother expecting an non.innov.SS TSmodel.")
 if  (is.null(filter$state) |  is.null(filter$track)) 
   {filter <- l(model,data, return.state=TRUE,return.track=TRUE)
    estimates <- filter$estimates
    filter    <- filter$filter
   }
 if (is.null(model$G))
  {m<-0
   G<-matrix(0,dim(model$F)[2],1)   # can't call compiled with 0 length arrays
   u <- matrix(0,nrow(filter$state),1)
  }
else
  {m <- dim(model$G)[2]
   G <-model$G
   u <- input.data(data)
  } 
sampleT  <-min(nrow(u), nrow(filter$state), dim(filter$track)[1])
 QQ <- model$Q %*% t(model$Q)           
 RR <- model$R %*% t(model$R) 
 n <- dim(model$F)[2]          
 if (compiled)
   {r<-.Fortran("smooth", 
                         state=filter$state,     # state, 
                         track=filter$track,     #trackerror
                         as.double(u),           # input
                         as.double(output.data(data)), # output
                         as.integer(n),                #n
                         as.integer(m),                #m
                         as.integer(dim(model$H)[1]),  #p 
                         sampleT =as.integer(sampleT), 
                         as.double(model$F),   
                         as.double(G),   
                         as.double(model$H), 
                         as.double(RR),
                         as.double(matrix(0,n,n)),   # scratch array
                         as.double(matrix(0,n,n)),   # scratch array
                         as.double(matrix(0,n,n)),   # scratch array
                         as.double(matrix(0,n,n)),   # scratch array
                         as.double(rep(0,n)),         # scratch array
                         DUP=.DSEDUP) [c("state","track")]
   }
 else   # S version
   {FF<-  model$F
    H <-  model$H
    #   filter$state is the one step ahead state estimate E{z(t)| y(t-1), u(t)}
    #   zt below is the filter state estimate E{z(t)| y(t), u(t+1)}
    #   filter$track is tracking error P(t|t-1)
    sm <- array(NA,dim(filter$state))  # smoother state estimate
    sm[sampleT,] <-filter$state[sampleT,]
    tr <- array(NA,dim(filter$track))  # smoother tracking error
    tr[sampleT,,] <-filter$track[sampleT,,]
    for (Time in (sampleT-1):1) 
      {K <- filter$track[Time,,] %*% t(H) %*% 
              solve(H %*% filter$track[Time,,] %*% t(H) + RR) #(A5)
       zt <- filter$state[Time,] + K %*% 
       		(output.data(data)[Time,] - H %*% filter$state[Time,]) #(A6)
       if (m!=0) zt <- zt - G %*% u[Time+1,]
       P <- filter$track[Time,,] - K %*% H %*% filter$track[Time,,]  #P(t|t)  (A7)
       P <- (P+t(P))/2 #force symmetry to avoid rounding problems
       J <- P %*% t(FF) %*% solve(filter$track[Time+1,,])             #(A8) 
       if (m==0)                                          
          sm[Time,] <- zt + J %*% (sm[Time+1,] - FF %*% zt)              #(A9) 
       else 
          sm[Time,] <- zt + J %*% (sm[Time+1,] - FF %*% zt - G %*% u[Time+1,]) #check   
       P <- P + J %*% (tr[Time+1,,] - filter$track[Time+1,,]) %*% t(J)    #(A10)
       tr[Time,,]<- (P+t(P))/2   #force symmetry to avoid rounding problems
      }
     r <- list(state=sm, track=tr) 
    }  # end S version
 
  state <- tframed(r$state, list(start=start(output.data(data)),
              frequency= frequency(output.data(data))), 
              names=dimnames(filter$state)[[2]]) 
  r <-(list(estimates=estimates, data=data, model=model, 
            filter=filter, smooth=list(state=state, track=r$track) ) )
   classed(r, "TSestModel") # constructor (smoother)
} # end of smoother
  



############################################################

#     parameter estimation functions   <<<<<<<<<<

############################################################


est.VARX <- function(...)
  {stop("est.VARX is defunct. Use est.VARX.ls or est.VARX.ar.")}

est.VARX.ls.old <- function(data, subtract.means=FALSE, standardize=FALSE,
    max.lag=NULL, trend=FALSE) 
{# Estimate VAR model with exogenous variable using lsfit(). 
 # Returns a TSestModel.
 # This is very similar to estimating Markov parameters (see est.SS.Mittnik).
 # Residuals,etc, are calculated by evaluating the estimated model with ARMA.
 # Data should be of class TSdata.
   if (is.null(max.lag)) max.lag <- 6
   data <- freeze(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   if(is.null(m))  m <- 0
   N <- periodsOutput(data)
   if (subtract.means)
    {if(m!=0)input.data(data)<-input.data(data)-t(matrix(apply(input.data(data),2, mean), m,N))
     output.data(data)<- output.data(data) - t(matrix(apply(output.data(data),2, mean), p,N))
    }
   if (standardize)
     {svd.cov <- svd(var(output.data(data)))
      output.data(data) <- output.data(data) %*% svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
     }
      # shift input to give feedthrough in one period
   if (m != 0) {z <- cbind(input.data(data)[2:N,],output.data(data)[1:(N-1),])}
   else z <- output.data(data)

 # The matrix Past is blocks of data:
 #  [in | out-1 | in-1 | out-2 | ... | in-max.lag | out-max.lag-1 ]
 # so the coef. matrix M has a corresponding structure.
 
   Past <- matrix(NA,N-max.lag,(p+m)*(max.lag))
   for (i in 0:(max.lag-1)) 
      Past[,(1+(m+p)*i):((m+p)*(1+i))] <-z[(max.lag-i):(N-1-i),]
   M <- t(lsfit(Past,output.data(data)[(max.lag+1):N,,drop=FALSE],intercept=trend)$coef)
   if (standardize && (m!=0))  # correct exogenous blocks for normalization
     {Tinv <- diag(svd.cov$d^.5, ncol=p)%*%svd.cov$u
      for (i in 0:(max.lag-1)) 
         M[,(1+(m+p)*i):(m+(m+p)*i)] <- Tinv %*% M[,(1+(m+p)*i):(m+(m+p)*i)]
     }
   TREND <- NULL
   if (trend)
     {TREND <- M[,1]
      M<-M[,2:(dim(M)[2])]
     }
   A <- array(NA, c(1 + max.lag, p, p))
   A[1,,] <-diag(1, p)
   for (i in 0:(max.lag-1)) 
          A[2+i,,] <-  -M[,(m+1+(m+p)*i):(m+p+(m+p)*i)]
   if (m==0) C <- NULL   # no exog. variable
   else               # NB. there is an implicit shift in input.data(data)
      {C <-array(NA,c(max.lag,p,m))
       for (i in 0:(max.lag-1)) 
         C[1+i,,] <- M[,(1+(m+p)*i):(m+(m+p)*i)]
      }
   B <- array(diag(1,p),c(1,p,p))
   l(ARMA(description="model estimated by est.VARX.ls",
        A=A,B=B,C=C,TREND=TREND, 
         input.names =  seriesNamesInput(data), 
        output.names = seriesNamesOutput(data)), data)
}


est.VARX.ls <- function(data, subtract.means=FALSE, re.add.means=TRUE,
   standardize=FALSE, unstandardize=TRUE, max.lag=NULL, 
   trend=FALSE, lag.weight=1.0, warn=TRUE) 
{# Estimate VAR model with exogenous variable using lsfit(). 
 # Returns a TSestModel.
 # This is very similar to estimating Markov parameters (see est.SS.Mittnik).
 # Residuals,etc, are calculated by evaluating the estimated model with ARMA.
 # Data should be of class TSdata.
 # lag.weight is an exponential weight applied to lags. It should be in (0,1].
   if (is.null(max.lag)) max.lag <- 6
   data <- freeze(data)
   names <- seriesNames(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   missing.data <- any(is.na(output.data(data)))
   if (0 != m) { # ensure same time window
	input.data(data) <- tfwindow(input.data(data),
              start.=startOutput(data), end.=endOutput(data), warn=FALSE)
	output.data(data) <- tfwindow(output.data(data),
              start.=startInput(data), end.=endInput(data), warn=FALSE)
	missing.data <- any(missing.data, is.na(input.data(data)))
   }
   N <- periodsOutput(data)
   if (standardize)
     {svd.cov <- svd(var(output.data(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if (subtract.means)
    {if(m!=0)
       {input.means<-apply(input.data(data),2, mean)
        input.data(data)<-input.data(data)-t(matrix( input.means, m,N))
       }
     output.means <- apply(output.data(data),2, mean)
     output.data(data)  <- output.data(data) - t(matrix(output.means, p,N))
    }
      # shift input to give feedthrough in one period
   if (m != 0) {z <- cbind(input.data(data)[2:N,],output.data(data)[1:(N-1),])}
   else z <- output.data(data)

 # The matrix Past is blocks of data:
 #  [in | out-1 | in-1 | out-2 | ... | in-max.lag | out-max.lag-1 ]
 # so the coef. matrix M has a corresponding structure.
 
   Past <- matrix(NA,N-max.lag,(p+m)*(max.lag))
   for (i in 0:(max.lag-1)) 
      Past[,(1+(m+p)*i):((m+p)*(1+i))] <-z[(max.lag-i):(N-1-i),] /(lag.weight^i)
   fit <- lsfit(Past,output.data(data)[(max.lag+1):N,,drop=FALSE],intercept=trend)
   if(missing.data)
     fit.res <-TSdata(output=fit$residual) #used only in the case of missing data for scaling
   M <- t(fit$coef)
   # correct exogenous blocks for normalization:
   if (standardize && (m!=0)) 
     {Tinv <- diag(svd.cov$d^.5, ncol=p)%*%svd.cov$u
      for (i in 0:(max.lag-1)) 
         M[,(1+(m+p)*i):(m+(m+p)*i)] <- Tinv %*% M[,(1+(m+p)*i):(m+(m+p)*i)]
     }
   TREND <- NULL
   if (trend)
     {TREND <- M[,1]
      M<-M[,2:(dim(M)[2]),drop=FALSE]
     }
   if (subtract.means & re.add.means)
    {if(m!=0) 
       {input.data(data)<-input.data(data) + t(matrix( input.means, m,N))
        Past<-Past +
           t(matrix(c(input.means,output.means),(p+m)*(max.lag),N-max.lag))
       }
     else
        Past <-Past+t(matrix(output.means, (p+m)*(max.lag),N-max.lag))
     output.data(data)  <- output.data(data) + t(matrix(output.means, p,N))
     # and correct estimation for non-zero mean:
     M<-M+est.VARX.mean.correction(Past, 
                       output.data(data)[(max.lag+1):N,,drop=FALSE],M, warn=warn)
    }
   A <- array(NA, c(1 + max.lag, p, p))
   A[1,,] <-diag(1, p)
   for (i in 0:(max.lag-1)) 
          A[2+i,,] <-  -M[,(m+1+(m+p)*i):(m+p+(m+p)*i),drop=FALSE] %*% 
                                 diag(lag.weight^i,p)
   if (m==0) C <- NULL   # no exog. variable
   else               # NB. there is an implicit shift in input.data(data)
      {C <-array(NA,c(max.lag,p,m))
       for (i in 0:(max.lag-1)) 
         C[1+i,,] <- M[,(1+(m+p)*i):(m+(m+p)*i),drop=FALSE]  %*% 
                                 diag(lag.weight^i,m)
      }
   B <- array(diag(1,p),c(1,p,p))
   model <-ARMA( description="model estimated by est.VARX.ls",
              A=A,B=B,C=C,TREND=TREND)
   seriesNames(model) <- seriesNames(data)

   if (standardize & unstandardize)
     {scalefac <- solve(scalefac)
      data <- scale(data,  scale=list(output=scalefac))
      model<- scale(model, scale=list(output=scalefac))
      if(missing.data) fit.res <- scale(fit.res, scale=list(output=scalefac))
     }

   if(missing.data)
       {if (warn) warning(
          "missing data kludge. Predictions are reconstructed from lsfit residuals")
        return(fake.TSestModel.missing.data(model,data, fit.res$output,max.lag))
       }
   else return(l(model, data, warn=warn))
}

fake.TSestModel.missing.data <- function(model,data, residual, max.lag)
  {pred <- rbind(matrix(NA,max.lag,nseriesOutput(data)),residual) + 
                                                           output.data(data)
   r <- list(estimates = residual.stats(pred, output.data(data), warn=warn),
               data = data, model = model)
   classed(r, "TSestModel") # fake missing data constructor
  }

est.VARX.mean.correction <- function(X, y, bbar,
                     fuzz=sqrt(.Machine$double.eps), warn=TRUE)
{# correction for model estimated with means subtracted
 Xbar <- t(array(apply(X, 2, mean), rev(dim(X))))
 ybar <- t(array(apply(y, 2, mean), rev(dim(y))))
 v <- svd(t(X)%*%X) # this is more robust than solve()
 if (warn && any(abs(v$d[1]*fuzz) > abs(v$d) ) ) 
   warning("The covariance matrix is nearly singular. Check for linearly related data.")
# if(1 == length(v$d))OmInv <- v$v %*% (1/v$d) %*% t(v$u)
# else OmInv <- v$v %*% diag(1/v$d) %*% t(v$u)	
#  following is equivalent
 OmInv <-  v$v %*% sweep(t(v$u),1,1/v$d, "*") 
 t(-OmInv %*% ( 
        (2*t(Xbar)%*%X - t(Xbar)%*%Xbar) %*% t(bbar) 
         - t(Xbar)%*%y - t(X)%*%ybar + t(Xbar)%*%ybar ))
}


est.VARX.ar <- function(data, subtract.means=FALSE,  re.add.means=TRUE,
      standardize=FALSE, unstandardize=TRUE, aic=TRUE, max.lag=NULL, 
      method="yule-walker", warn=TRUE) 
{
# Estimate VAR model with exogenous variable using ar(). Returns a TSestModel.
# Residuals,etc, are calculated by evaluating the estimated model with ARMA.
# Use Splus procedure ar  and combine exogeneous variables.
# Note: ar uses a Yule-Walker approach (uses autocorrelations) so effectively the 
#   model is for data with means removed. Thus subtract.means does not make much
#   difference and re.add.means must be T to get back to a model for the 
#   original data.
# Conventon for AR(0)  and sign are changed to ARMA format.
# Data should be of class TSdata.
# The exog. variable is shifted so contemporaneous effects enter.
# the model for the exog. variable (as estimated by ar() is  discarded.
   data <- freeze(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   N <- periodsOutput(data)
   if (standardize)
     {svd.cov <- svd(var(output.data(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if(m!=0) input.means<-apply(input.data(data),2, mean)
   output.means <- apply(output.data(data),2, mean)
   if (subtract.means)
    {if(m!=0) input.data(data)<-input.data(data)-t(matrix( input.means, m,N))
     output.data(data)  <- output.data(data) - t(matrix(output.means, p,N))
    }
   if (m==0)  zdata <- output.data(data)  # no exog. variable
   else       zdata <- cbind(input.data(data),output.data(data))
         # NB. there is an implicit shift in input.data(data) in the line above
         #   because ar estimates lag parameters,
         # so C[1,,] is for one lag in u.
   if (is.null(max.lag))  AC<- DSE.ar(zdata, aic=aic, method=method)
   else                   AC<- DSE.ar(zdata, aic=aic, method=method, order.max=max.lag)
   if (re.add.means)
     {if (subtract.means)
        {if(m!=0) input.data(data)<-input.data(data) + t(matrix( input.means, m,N))
         output.data(data)  <- output.data(data) + t(matrix(output.means, p,N))
        }
      # and correct estimation for non-zero mean:
      max.lag <- AC$order
      if (max.lag==0) stop("all lags eliminated by AIC order selection.")
      if (m==0)  zdata <- output.data(data)  # no exog. variable
      else       zdata <- cbind(input.data(data),output.data(data))
      Past <- matrix(NA,N-max.lag,(p+m)*(max.lag))
      for (i in 0:(max.lag-1)) 
         Past[,(1+(m+p)*i):((m+p)*(1+i))] <-zdata[(max.lag-i):(N-1-i),]
      M <-est.VARX.mean.correction(Past, zdata[(max.lag+1):N,,drop=FALSE],
                  matrix(aperm(AC$ar, c(2,3,1)),m+p,(m+p)*max.lag), warn=warn)
      AC$ar<-AC$ar + aperm(array(M,c(m+p,m+p,max.lag)), c(3,1,2))
     }
   A <- array(0, c(1 + AC$order, p, p))
   A[1,,] <-diag(1, p)
   if (0==AC$order)
      warning("lagged output variables eliminated by AIC order selection.")
   else 
      A[2:(1+AC$order),,] <-  -AC$ar[, (m+1):(m+p), (m+1):(m+p), drop=FALSE]
   if (m==0) C <- NULL
   else 
     {if (0==AC$order)
             {warning("input variables eliminated by AIC order selection.")
              C <-array(0, c(1,p,m))
             }
      else C <-array(AC$ar[,(m+1):(m+p),1:m],c(AC$order,p,m))
     }
   B <- array(diag(1,p),c(1,p,p))
   model <-ARMA(
            description="model estimated by est.VARX.ar",
            A=A,B=B,C=C, 
            names = seriesNames(data) )
 
   if (standardize & unstandardize)
     {scalefac <- solve(scalefac)
      data <- scale(data,  scale=list(output=scalefac))
      model<- scale(model, scale=list(output=scalefac))
     }
   l(model, data, warn=warn)
}

old.est.VARX.ar <- function(data, subtract.means=FALSE,  re.add.means=TRUE,
    standardize=FALSE, unstandardize=TRUE, aic=TRUE, max.lag=NULL, 
    method="yule-walker") 
{
# Estimate VAR model with exogenous variable using ar(). Returns a TSestModel.
# Residuals,etc, are calculated by evaluating the estimated model with ARMA.
# Use Splus procedure ar  and combine exogeneous variables.
# Conventon for AR(0)  and sign are changed to ARMA format.
# Data should be of class TSdata.
# The exog. variable is shifted so contemporaneous effects enter.
# the model for the exog. variable (as estimated by ar() is  discarded.
   data <- freeze(data)
   m <-  nseriesInput(data)
   p <- nseriesOutput(data)
   if(is.null(m))  m <- 0
   N <- periodsOutput(data)
   if (subtract.means)
    {if(m!=0)
       {input.means<-apply(input.data(data),2, mean)
        input.data(data)<-input.data(data)-t(matrix( input.means, m,N))
       }
     output.means <- apply(output.data(data),2, mean)
     output.data(data)  <- output.data(data) - t(matrix(output.means, p,N))
    }
   if (standardize)
     {svd.cov <- svd(var(output.data(data)))
      scalefac <- svd.cov$u %*% diag(1/svd.cov$d^.5, ncol=p)
      data <- scale(data, scale=list(output=scalefac))
     }
   if (m==0)    # no exog. variable
      {if (is.null(max.lag))   AC<- DSE.ar(output.data(data), aic=aic, method=method)
       else AC<- DSE.ar(output.data(data), order.max=max.lag, aic=aic, method=method)
       A <- array(0, c(1 + AC$order, p, p))
       A[1,,] <-diag(1, p)
       if (0<AC$order) A[2:(1+AC$order),,] <-  -AC$ar
       C <- NULL
      }
   else               # NB. there is an implicit shift in input.data(data)
      {               #   because ar estimates lag parameters,
                      # so C[1,,] is for one lag in u.
       if (is.null(max.lag))
         AC<- DSE.ar(cbind(input.data(data),output.data(data)), aic=aic, method=method) 
       else   
         AC<- DSE.ar(cbind(input.data(data),output.data(data)), order.max=max.lag, aic=aic, method=method) 
       A <- array(0, c(1 + AC$order, p, p))
       A[1,,] <-diag(1, p)
       if (0<AC$order) 
         {A[2:(1+AC$order),,] <-  -AC$ar[, (m+1):(m+p), (m+1):(m+p), drop=FALSE]
          C <-array(AC$ar[,(m+1):(m+p),1:m],c(AC$order,p,m))
         }
      }
#   B <- array(AC$var.pred,c(1,p,p))
   B <- array(diag(1,p),c(1,p,p))
   model <-ARMA(description="model estimated by est.VARX.ls",
                   A=A,B=B,C=C, 
                  names = seriesNames(data) )

   if (standardize & unstandardize)
     {scalefac <- solve(scalefac)
      data <- scale(data,  scale=list(output=scalefac))
      model<- scale(model, scale=list(output=scalefac))
     }
   if (subtract.means & re.add.means)
    {if(m!=0) input.data(data)<-input.data(data) + t(matrix( input.means, m,N))
     output.data(data)  <- output.data(data) + t(matrix(output.means, p,N))
    }
   l(model, data)
}

est.SS.from.VARX <- function(data, warn=TRUE, ...){
	# estimate a VARX model and convert to state space
	#  estimate a nested-balanced state space model by svd a la Mittnik from
	#   least squares estimate of VAR coefficents.
	model <-est.VARX.ls(data, warn=warn, ...)
	l(to.SS(model$model), model$data, warn=warn)
}




############################################################

#     model balancing functions   <<<<<<<<<<

############################################################

balance.Mittnik <- function(model, n=NULL){
 # calculate a nested-balanced state space model by svd
   # n is intended primarily for producing a state space model from the markov
#  parameters of an ARMA model, but if it is supplied with an SS model the
#  result will be a model with state dimension n based on the n largest
#  singular values of the svd of a Hankel matrix of markov parameters generated
#  by the original model. If n is not supplied then the singular values are
#  printed and the program prompts for n.

  if (is.TSestModel(model)) model <- TSmodel(model)
  if  (!is.TSmodel(model)) stop("balance.Mittnik expecting a TSmodel.")
  m <- nseriesInput(model)
  newmodel <- balance.Mittnik.svd(markov.parms(model), m, n)$model
  newmodel$description <- paste(model$description,"converted to", newmodel$description)
  seriesNames(newmodel)  <-  seriesNames(model)
  newmodel
}


balance.Mittnik.svd <- function(M, m, n=NULL)
{ # Calculate a nested-balanced state space model by svd a la Mittnik.
  # If state dim n is supplied then svd criteria are not calculated and 
  # the given n is used. Otherwise, the singular values are printed and 
  # the program prompts for n.
  # M is a matrix with p x (m+p)  blocks giving the markov parameters,
  # that is, the first row of the Hankel matrix. It can be generated from the
  # model as in the function markov.parms, or from the data, as in the function
  # est.SS.Mittnik.
  # m is the dimension of input series, which is needed to decompose M.
  # The output dimension p is taken from nrow(M).

#                   # Form k block Hankel Matrix from M.
   p <- nrow(M)     # dim of endo. series
   r <- m + p       # each sub-matrix is p x r   (= p x (m+p) )
   k<- dim(M)[2] / r
   Hkl <- matrix(0, p*k, r*k) 
   for(i in 1:k)             # Hankel with 0s in SE half  
     Hkl[(1+p*(i-1)):(p*i), 1:(r*(1+k-i))]<- M[,(1+r*(i-1)):(r*k),drop=FALSE]
# note: if last block of M is 0, which is often the case, then filling
#   with zeros, as above, is correct.
#   k <- k %/% 2
#   Hkl <- matrix(0, p*k, r*k) 
#   for (i in 1:k)              # (smaller) completely filled Hankel
#     {for (j in 1:k)       
#        {Hkl[(1+p*(i-1)):(p*i), (1+r*(j-1)):(r*j)] <-M[ ,(1+r*(i+j-2)):(r*(i+j-1))] }}

   svd.of.Hkl <- svd(Hkl)
   shifted.Hkl <- Hkl[,(r+1):dim(Hkl)[2]]
   shifted.Hkl <- cbind(shifted.Hkl,matrix(0,p*k,r))
   rtQ <- diag(sqrt(svd.of.Hkl$d))
   rtQ.inv <- diag(1/sqrt(svd.of.Hkl$d))   
   o.mtr <- svd.of.Hkl$u %*% rtQ 
   o.inv <- t(svd.of.Hkl$u %*% rtQ.inv )
   c.mtr <- rtQ %*% t(svd.of.Hkl$v)
   c.inv <- t(rtQ.inv %*% t(svd.of.Hkl$v))
   svd.crit <- svd.of.Hkl$d
   c.mtr[svd.crit==0,] <-0     # this gets rid of NAs but the resulting
   c.inv[,svd.crit==0] <-0     #  model may not be minimal
   o.inv[svd.crit==0,] <-0     # eg. F may have zero rows and columns
   o.mtr[,svd.crit==0] <-0     #   (specifically if the model included a trend)
   if (is.null(n))
     {svd.crit <- svd.criteria(svd.of.Hkl$d)
      n <- read.int("Enter the number of singular values to use for balanced model: ")
     }
   H  <- o.mtr[1:p,1:n,drop=FALSE]
   FF <- o.inv[1:n,,drop=FALSE] %*% shifted.Hkl%*%c.inv[,1:n,drop=FALSE]
   if (m == 0) G <- NULL             # no exog.
   else        G <- c.mtr[1:n,1:m,drop=FALSE]  # exog.
   K  <- c.mtr[1:n,(m+1):(m+p),drop=FALSE]
   FF <- FF + K%*%H  #converts from system with lagged outputs as 
                    #   inputs see Mittnik p1190.
  #browser()
  #  good checks here are (with n=length(svd.of.Hkl$d)):
  #   0 == max(abs(shifted.Hkl - o.mtr[,1:n] %*% (FF-K%*%H) %*% c.mtr[1:n,]))
  #   0 != max(abs(shifted.Hkl))
  model<-SS(description="nested-balanced model a la Mittnik",
              F=FF,G=G,H=H,K= K)
  list(crit=svd.crit,model=model)
}

############################################################

#     statistical test functions  <<<<<<<<<<

############################################################

Portmanteau <- function(res){
  # Portmanteau statistic for residual
  if (is.R()) if (!require("ts", warn.conflicts = FALSE)) stop("package ts is required.")
  ac <- acf(res,type="covariance", plot=FALSE)$acf
  p <- dim(ac)[1]
#  a0 <- solve(ac[1,,])  the following is more robust than solve for
#              degenerate densities
  v <- svd(ac[1,,])
#  if(1 == length(v$d)) a0 <- v$v %*%     (1/v$d) %*% t(v$u)
#	           else a0 <- v$v %*% diag(1/v$d) %*% t(v$u)	
#  following is equivalent
  a0 <-  v$v %*% sweep(t(v$u),1,1/v$d, "*") 
  P <-0
  for (i in 2:p) 
    {P <- P + sum(diag(t(ac[i,,]) %*% a0 %*% ac[i,,] %*% a0))
    }
  dim(res)[1]*P
}


check.residuals <- function(obj, ...)  UseMethod("check.residuals")
# autocorrelations <- function(obj, ...) UseMethod("check.residuals")

check.residuals.TSestModel <- function(obj, ...){
   invisible(check.residuals(obj$estimates$pred - output.data(obj), ...))}

check.residuals.TSdata <- function(obj, ...){
	invisible(check.residuals(output.data(obj), ...)) }

check.residuals.default <- function(obj, ac=TRUE, pac=TRUE,
     select=seq(nseries(obj)), drop=NULL, plot.=TRUE, 
     graphs.per.page=5, verbose=FALSE)
{#  select can indicate column indices of  residuals to use.
 # If drop is not null it should be a vector of the row indices of residuals
 #  which are not to be used. Typically this can be used to get rid
 #  of bad initial conditions (eg. drop=seq(10) ) or outliers.
  resid <- select.series(obj, series=select)
  if (!is.null(drop))  resid <- resid[-drop,, drop=FALSE]
  mn<-apply(resid,2,mean)
  acr <-NULL
  pacr<-NULL
  cov <- var(resid)
  if (verbose) 
    {cat("residual means: \n") ; print(mn); cat("\n")
     cat("residual cov  : \n") ; print(cov) ; cat("\n")
    }
  p <- nseries(resid)
  resid0 <- sweep(resid, 2, mn, FUN="-")
#  resid0 <- resid - t(array(apply(resid,2,mean),rev(dim(resid)))) # mean 0
  cusum <- apply(resid0,2,cumsum)/ t(array(diag(var(resid0)),rev(dim(resid0))))
  if (is.R()) if (!require("ts", warn.conflicts = FALSE)) stop("package ts is required.")
  if(plot. &&  dev.cur() != 1 ) 
    {graphs.per.page <- min(p, graphs.per.page)
     names <- seriesNames(resid)
     old.par <-par(mfcol = c(3, graphs.per.page), mar = c(5.1, 4.1,3.1, 0.1) ) #c(5,4.1,5,0.1) c(2.1, 4.1,3.1, 0.1)
     on.exit(par(old.par))
     for (i in 1:p) 
          {tfOnePlot(resid[,i], xlab=names[i])
           if (i %% graphs.per.page == ceiling(graphs.per.page/2)) 
                 title(main ="Residuals ")
           if (i %% graphs.per.page == 0) 
                 title(main = paste("       page ", floor(i / graphs.per.page)))
           tfOnePlot(cusum[,i], xlab=names[i])
           if (i %% graphs.per.page == ceiling(graphs.per.page/2)) 
                 title(main = "Cusum")
           if (exists("ksmooth")) 
                       rd<-ksmooth(resid[,i],bandwidth=var(resid[,i])^0.5)
           else  if (exists("density"))
                       rd <- density(resid[,i],       bw=var(resid[,i])^0.5)
           else
        stop("Neither ksmooth nor density are available to calculate the plot.")
           plot(rd,type="l",xlab=names[i],ylab="")
           if (i %% graphs.per.page == ceiling(graphs.per.page/2))
              title(main = "kernel estimate of residual distribution")
          }
     if (ac)
       {par(mfrow = c(1, 1), mar = c(2.1, 4.1,3.1, 0.1), oma=c(0,0,5,0) )
        acr  <-acf(resid, plot=TRUE)$acf
        mtext("Autocorrelations", side=3, outer=TRUE, cex=1.5)
       }
     if (pac)
       {par(mfrow = c(1, 1), mar = c(2.1, 4.1,3.1, 0.1), oma=c(0,0,5,0) )
        pacr <-acf(resid, plot=TRUE, type= "partial")$acf 
        mtext("Partial Autocorrelations", side=3, outer=TRUE, cex=1.5)
       }
    }
  else 
    {if (ac)  acr  <-acf(resid, plot=FALSE)$acf
     if (pac) pacr <-acf(resid, plot=FALSE, type= "partial")$acf 
    }
  if (ac & verbose)
    {cat("residual auto-correlations:\n")
     cat("      lag:   ")
     for (i in 1: dim(acr)[1])  cat(i,"        ")  
     cat("\n") 
     for (i in 1:p) { cat(i,": "); cat(acr[,i,i]); cat("\n")  }
    }
  if (pac & verbose)
    {cat("partial auto-correlations:\n")
     cat("      lag:   ")
     for (i in 1: dim(pacr)[1])  cat(i,"        ")  
     cat("\n") 
     for (i in 1:p) { cat(i,": "); cat(pacr[,i,i]); cat("\n")  }
    }
#  cat("residual normality tests:\n")
#  cat("hetros. tests:\n")

  skewness <- apply(sweep(resid,2,mn)^3,2,mean) / diag(cov)^(3/2)
  kurtosis <- apply(sweep(resid,2,mn)^4,2,mean) / diag(cov)^2
  invisible(list(residuals=resid, mean=mn, cov=cov, acf=acr, pacf=pacr, 
                 cusum=cusum, skewness=skewness, kurtosis=kurtosis))
}

information.tests <- function(..., sample.start=1,sample.end=NULL,
		 Print=TRUE, warn=TRUE){
 # print model selection criteria
     #  for models statistics ..., where ... are the names of the
  #  list information as returned by like.
  #  likes returns neg. log likelihood as lst$like[4].
  if (Print) criteria.table.heading()
  values <- NULL
  options(width=100)
  for (lst in list(...) ) 
    {z <-information.tests.calculations(lst,
    		 sample.start=sample.start,sample.end=sample.end, warn=warn)
     values <-rbind(values, z)
     if (Print)
       {print(c(z),digits=4)
        cat("\n")
       }
    }
  if (Print & (1 <dim(values)[1]) )
    {cat("opt     ")
     opt <-apply(values,2,order)[1,]  # minimum
     for (i in 1:length(opt)) cat(opt[i],"      ")
     cat("\n")
    }
  if (Print) criteria.table.legend()
  invisible(values)
}


information.tests.calculations <- function # return model selection criteria
     (lst, sample.start=1,sample.end=NULL, warn=TRUE){
   resid <- lst$estimates$pred-output.data(lst$data)
    # the following line is just to work around a bug with old style time series
   if (ncol(output.data(lst$data))==1) dim(resid) <- dim(output.data(lst$data))
   if (is.null(sample.end)) sample.end <- nrow(resid)
   resid <- resid[sample.start:sample.end,,drop=FALSE]
   ml <-   residual.stats(resid, NULL, warn=warn)$like[1] # neg.log(likelihood)
# previously   ml <-   L(resid)[1]     # neg. log( likelihood ).
   n  <-  length(coef(lst$model))      # No. of parameters.
   # nt is theorical dimension of parameter space    n(m+2p)
   if (is.ARMA(lst$model)) nt <- NA 
   if (is.SS(lst$model))  
      if (is.null(lst$model$G)) nt <- nrow(lst$model$F)*2*nrow(lst$model$H)
      else   nt <- nrow(lst$model$F)*(ncol(lst$model$G)+2*nrow(lst$model$H))
   r  <- nrow(lst$estimates$pred)*ncol(lst$estimates$pred) #No. of residuals.
   port <-Portmanteau(resid)
   aic <-  2*ml + 2*n             # AIC
   ops  <- options(warn=-1)
   on.exit(options(ops))
   bic <-  2*ml + n*log(r)        # BIC
   gvc <-  2*ml - 2*r*log(1-n/r)  # GCV
   rice<-  2*ml - r*log(1-2*n/r)  # RICE
   fpe <-  2*ml + r*(log(1+(n/r))-log(1-(n/r)))  # FPE
   taic <-  2*ml + 2*nt             # AIC
   tbic <-  2*ml + nt*log(r)        # BIC
   tgvc <-  2*ml - 2*r*log(1-nt/r)  # GCV
   trice<-  2*ml - r*log(1-2*nt/r)  # RICE
   tfpe <-  2*ml + r*(log(1+(nt/r))-log(1-(nt/r)))  # FPE
   z<- matrix(c(port,ml,aic,bic,gvc,rice,fpe,taic,tbic,tgvc,trice,tfpe),1,12)
   dimnames(z)<-list(NULL,c("port","like","aic","bic","gvc","rice","fpe","taic",
                   "tbic","tgvc","trice","tfpe"))
   z
}

criteria.table.heading <- function(){
        cat("                    based on no.of parameters     based on theoretical parameter dim.","\n")
        cat("       PORT  -ln(L)  AIC   BIC   GVC  RICE   FPE   AIC   BIC   GVC   RICE   FPE\n")
}

criteria.table.nheading <- function(){
        cat("                    based on no.of parameters     based on theoretical parameter dim.","\n")
        cat("dim.   PORT  -ln(L)  AIC   BIC   GVC  RICE   FPE   AIC   BIC   GVC   RICE   FPE\n")
}

criteria.table.legend <- function(){
        cat("  PORT  - Portmanteau test                     ") 
        cat("  -ln(L)- neg. log likelihood                  \n")
        cat("  AIC   - neg. Akaike Information Criterion    ")  
        cat("  BIC   - neg. Bayes  Information Criterion    \n")
        cat("  GVC   - Generalized Cross Validation         ") 
        cat("  RICE  - Rice Criterion                       \n")
        cat("  FPE   - Final Prediction Error               \n") 
        cat(" WARNING - These calculations do not account for trend parameters in ARMA models.\n") 
}


svd.criteria <- function(sv){
   cat("\n SINGULAR VALUES OF THE HANKEL MATRIX:\n")
   print(sv,digits=3)
   svsqr.vct <- sv^2/sum(sv^2)    
   cat("\n GELFAND & YAGLOM INFORMATION CRITERIA:\n")
   GY <-cumsum(log(1-svsqr.vct))/sum(log(1-svsqr.vct))
   print(GY,digits=3)
   cat("\n RATIO OF SINGULAR VALUES TO MAX SINGULAR VALUE:\n")
   print(sv/sv[1],digits=3)
   cbind(sv,GY,sv/sv[1])
}   

#######################################################################

tfwindow.TSdata <- function(x, start=NULL, end=NULL, tf=NULL, warn=TRUE)
{# window a TSdata object
  if (0 != nseriesInput(x))  input.data(x) <- 
      tfwindow(input.data(x), start=start, end=end, tf=tf, warn=warn)
  if (0 != nseriesOutput(x)) output.data(x) <- 
      tfwindow(output.data(x), start=start, end=end, tf=tf, warn=warn)
  x
}

window.TSdata <- tfwindow.TSdata


combine <- function(e1,e2)UseMethod("combine")
combine.default <- function(e1,e2){list(e1,e2)}


combine.TSdata <- function(e1,e2)
{# make a new TSdata object with the two objects
 output.data(e1) <- tbind(output.data(e1),output.data(e2)) 
 if ((0 != (nseriesInput(e1))) & (0 != (nseriesInput(e2))))
       input.data(e1) <- tbind(input.data(e1),input.data(e2)) 
 else  {if (0 != (nseriesInput(e2)))  input.data(e1) <-input.data(e2)  }
 e1
}


trim.na.TSdata <- function(x, start.=TRUE, end.=TRUE)
{# trim NAs from the ends of TSdata.
 # (Observations for all series are dropped if any one contains an NA.)
 # if start.=F then beginning NAs are not trimmed.
 # If end.=F   then ending NAs are not trimmed.
 # The same truncation is applied to both input and output
 p <- nseriesOutput(x)
 m <- nseriesInput(x)
 if (m==0)
   mat <- trim.na(output.data(x))
 else
   mat <- trim.na(tbind(input.data(x),output.data(x)),start.=start.,end.=end.)
 tf <- tframe(mat)
 if (m!=0)
   {sn <- seriesNamesInput(x)
    input.data(x)  <- tframed(mat[,1:m, drop=FALSE], tf) 
    seriesNamesInput(x) <- sn
   }
 sn <- seriesNamesOutput(x)
 output.data(x) <- tframed(mat[,(m+1):(m+p), drop=FALSE], tf)
 seriesNamesOutput(x) <- sn
 x
}



diff.log <- function(x,  lag = 1, base = exp(1))
{#Calculate the difference from lag periods prior for log of data.
diff(log(x, base = base), lag = lag)
}


ytoypc <- function(ser) {
  # Convert level data to year over year percent change.
  # note: percent.change can alter the name, so grab it first.
  nm <- paste("y to y %ch", seriesNames(ser))
  ser <- percent.change(ser, lag=frequency(ser))
  seriesNames(ser) <- nm
  ser
 }
 


percent.change <- function(obj, ...) UseMethod("percent.change")

percent.change.list <- function(obj, ..., base=NULL, lag=1, cumulate=FALSE, e=FALSE)
  {#Calculate the percent change relative to the data lag periods prior.
   #... should be a list of objects to which percent.change can be applied.
   pchange <- list(percent.change(obj, base=base, lag=lag, cumulate=cumulate, e=e))
   for (mat in list(...))
          pchange <- append(pchange,list(
	    percent.change(mat, base=base, lag=lag, cumulate=cumulate, e=e)))
   pchange
  }

percent.change.default <- function(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE)
{#Calculate the percent change relative to the data lag periods prior.
 # obj should be  a  matrix or vector.
 # If cumulate is TRUE then the data is cumulated first. cumulate can be
 # a logical vector with elements corresponding to columns of m.
 # If e is T the exponent of the series is used (after cumulating 
 #   if cumulate is T). e can be
 # a logical vector with elements corresponding to columns of m.
 # If base is provided it is treated as the first period value 
 #  (prior to differencing). It is prefixed to the m prior to 
 #  cumulating. It should be a vector of length dim(m)[2]. 
 #  (If e is TRUE then base should be log of the original data).
   cls <- dseclass(obj)
   # note next has to be applied to a shorter object in the end
   if (is.tframed(obj)) tf <- list(end=end(obj), frequency=frequency(obj))
   else tf <- NULL
   if (is.null(dim(obj)))
     {vec <- TRUE
      obj <- matrix(obj, length(obj),1)
     }
   else vec <- FALSE
   mm <- rbind(base,obj)
   if (any(cumulate))
          mm[,cumulate] <-apply(mm[,cumulate,drop=FALSE],2,cumsum)
   if (any(e)) mm[,e] <- exp(mm[,e,drop=FALSE])
   N <- dim(mm)[1]
   pchange <-100*(mm[(lag+1):N,,drop=FALSE] - 
                    mm[1:(N-lag),,drop=FALSE])/mm[1:(N-lag),,drop=FALSE]
   if (vec) pchange <- pchange[,1]
   dseclass(pchange) <- cls
   if (!is.null(tf)) tframed(pchange, tf) else pchange
}

percent.change.TSestModel <- function(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE)
  {#The percent change calculation is done
   # with $estimates$pred and the result is an object of class TSdata
   TSdata(output=percent.change(obj$estimates$pred))
  }

percent.change.TSdata <- function(obj, base=NULL, lag=1, cumulate=FALSE, e=FALSE)
  {# The percent change calculation is done
   # with input and output and the result is an object of class TSdata.
   if (0 != (nseriesInput(obj)))  input.data(obj)  <- percent.change(input.data(obj))
   if (0 != (nseriesOutput(obj))) output.data(obj) <- percent.change(output.data(obj))
   obj
  }


standardize <- function(ser){
	if (!is.matrix(ser)) stop("series should be a matrix.")
	means <- apply(ser, 2, mean)
	new <- ser - t(matrix(means, ncol(ser), nrow(ser)))
	svd.cov <- svd(var(new))
        scalefac <- svd.cov$u %*% diag(1/svd.cov$d^0.5, ncol = length(svd.cov$d))
        new <- new %*% t(scalefac)
	tframed(new, tf=tframe(ser), names=seriesNames(ser))
    }


############################################################

#     model and data scaling functions   <<<<<<<<<<

############################################################

# the following makes a generic copy of the function in the library to resolve
#  problems with version of S which do not have a generic function.
#  ( see also .First.lib at the beginning of dse1a.s
 # otherwise the following produces warning messages
invisible(
   if (!exists("scale.default"))
     {if (exists("scale")) scale.default <- scale  
      scale <- function(x, center=FALSE, scale = TRUE) UseMethod("scale")
     }
  )
  


scale.TSdata <- function(x, center=FALSE, scale=NULL) 
{# scale should be a list with two matrices or vectors, named input and output,
 # giving the multiplication factor for inputs and outputs.
 # vectors are treated as diagonal matrices.
 # If input or output are NULL then no transformation is applied.
 # The resulting data has inputs and outputs which are different from
 #  the original in proportion to scale. ie. if S and T are output and input
 #  scaling matrices then 
 #          y'(t) = S y(t) where y' is the new output
 #          u'(t) = S u(t) where u' is the new input
 sc <- scale$input # input.data(scale) causes problems if scale is a vector
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("input scale must be a vector or matrix")
    d <- input.data(x)
    tf <- tframe(d)
    names <- seriesNames(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    seriesNames(d) <- names
    input.data(x) <- d 
   }
 sc <- scale$output # output.data(scale) causes problems if scale is a vector
 if (!is.null(sc))
   {if (! (is.matrix(sc) | is.vector(sc))) stop("output scale must be a vector or matrix")
    d <- output.data(x)
    tf <- tframe(d)
    names <- seriesNames(d)
    if(is.matrix(sc))      d <- d %*% t(sc)
    else if(1==length(sc)) d <- d * sc 
    else                   d <- d %*% diag(sc)                             
    tframe(d) <- tf
    seriesNames(d) <- names
    output.data(x) <- d 
   }
 x
}


scale.TSestModel <- function(x, center=FALSE, scale=NULL)
   {scale(TSmodel(x), center=center, scale=scale)}

checkScale <- function(x, scale) UseMethod("checkScale")

checkScale.TSestModel <- function(x, scale){checkScale(TSmodel(x), scale)}

checkScale.TSmodel <- function(x, scale) 
{# This function only checks for some error conditions.
 # input.data(scale) and output.data(scale) cause problems if scale is a vector
 if (!is.null(scale$input))
   {if (is.matrix(scale$input))
       {if (any(svd(scale$input)$d == 0))  
        stop("scale input transformations must be non singular.")
       }
    else if(any(scale$input == 0)) stop("input scale elements must be non zero.")
   }
 if (!is.null(scale$output))
   {if (is.matrix(scale$output))
      {if (any(svd(scale$output)$d == 0))  
       stop("output scale transformations must be non singular.")
      }
    else if(any(scale$output == 0)) stop("output scale elements must be non zero.")
   }
 invisible(TRUE)
}


scale.innov <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 if (!is.null(scale$input)) 
   {sc <- scale$input
    if(is.vector(sc)) sc <- diag(sc, nseriesInput(x))
    x$G <- x$G %*% solve(sc)
   }
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))          
 x$H <- sc %*% x$H
 x$K <- x$K %*% solve(sc)
 set.parameters(x)
}

scale.non.innov <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 if (!is.null(scale$input)) 
   {sc <- scale$input
    if(is.vector(sc)) sc <- diag(sc, nseriesInput(x))
    x$G <- x$G %*% solve(sc)
   }
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))          
 x$H <- sc %*% x$H
 x$R <- sc %*% x$R %*% solve(sc)
 set.parameters(x)
}

scale.ARMA <- function(x, center=FALSE, scale=NULL)
{if (!checkScale(x, scale)) stop("scaling error.")
 sc <- scale$output
 if(is.vector(sc)) sc <- diag(sc, nseriesOutput(x))   
 x$A <- polyprod(sc, polyprod(x$A, solve(sc)))
 x$B <- polyprod(sc, polyprod(x$B, solve(sc)))
 if (!is.null(x$C)) 
   {sci <- scale$input
    if (!is.null(sci)) 
       {if(is.vector(sci)) sci <- diag(sci, nseriesInput(x))          
        x$C <- polyprod(x$C, solve(sci))
       }
    x$C <- polyprod(sc, x$C)
   }
 if (!is.null(x$TREND))  x$TREND <- sc %*% x$TREND
 set.parameters(x)
}


############################################################

#     Methods which are generic for models and TSdata    <<<<<<<<<<

# generic methods for TS structures (ie with input and output) <<<<<<<<<<
#                   see also tframe.s                          <<<<<<<<<<

############################################################################
# periods and periods.default are defined in tfame.s

periodsInput <- function(x)UseMethod( "periodsInput")
periodsOutput <- function(x)UseMethod("periodsOutput")

startInput <- function(x)UseMethod("startInput")
startOutput <- function(x)UseMethod("startOutput")

endInput <- function(x)UseMethod("endInput")
endOutput <- function(x)UseMethod("endOutput")

frequencyInput <- function(x)UseMethod("frequencyInput")
frequencyOutput <- function(x)UseMethod("frequencyOutput")




input.data <- function(x,  series=seqN(nseriesInput(x)))UseMethod("input.data")
output.data <- function(x, series=seqN(nseriesOutput(x)))UseMethod("output.data")

input.data.default <- function(x, series=seqN(nseriesInput(x)))
   select.series(x$input, series=series)
output.data.default <- function(x, series=seqN(nseriesOutput(x)))
   select.series(x$output, series=series)

"input.data<-" <- function(x,  value)  UseMethod("input.data<-")
"output.data<-" <- function(x,  value) UseMethod("output.data<-")

"input.data<-.default" <- function(x,  value) {x$input <- value; x}
"output.data<-.default" <- function(x,  value){x$output <- value; x}

# The logic (revised as of May26, 1998) is that seriesNames can be an attribute
# of any object (like a matrix) and has been moved to tframe.



 seriesNamesInput <- function(x)UseMethod( "seriesNamesInput")
seriesNamesOutput <- function(x)UseMethod("seriesNamesOutput")

 "seriesNamesInput<-" <- function(x, value)UseMethod( "seriesNamesInput<-")
"seriesNamesOutput<-" <- function(x, value)UseMethod("seriesNamesOutput<-")


############################################################################

#    methods for TSmodels  <<<<<<<<<<
# and check elsewhere too

############################################################################

seriesNames.TSmodel <- function(x)
  {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

seriesNamesOutput.TSmodel <- function(x)
 {# return output names if available in the object,
  # otherwise return "out" pasted with integers.
  if (!is.null(attr(x, "seriesNamesOutput")))
                          return(attr(x, "seriesNamesOutput"))
  else if (0 != nseriesOutput(x)) 
                   return(paste("out", seq(nseriesOutput(x)), sep=""))  
  else return(NULL)
 }

seriesNamesInput.TSmodel <- function(x)
 {# return input names if available in the object,
  # otherwise return "in" pasted with integers.
  if (!is.null(attr(x, "seriesNamesInput")))
                          return(attr(x, "seriesNamesInput"))
  else if (0 != nseriesInput(x)) 
                   return(paste("in", seq(nseriesInput(x)), sep=""))  
  else return(NULL)
 }

"seriesNames<-.TSmodel" <- function(x, value)
   { seriesNamesInput(x) <-  value$input
    seriesNamesOutput(x) <-  value$output
    x
   }

"seriesNamesInput<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != nseriesInput(x))
       stop("model input dimension and number of names do not match.")
    attr(x, "seriesNamesInput")  <- value
    x
   }

"seriesNamesOutput<-.TSmodel" <- function(x,  value)
   {if (!is.null(value) && length(value) != nseriesOutput(x))
       stop("model output dimension and number of names do not match.")
    attr(x, "seriesNamesOutput")  <- value
    x
   }


############################################################################

#    methods for TSestModels  <<<<<<<<<<
# and check elsewhere too   (esp. for start, end and frequency)

############################################################################
start.TSestModel <- function(x)start(x$data)
startInput.TSestModel <- function(x)startInput(x$data)
startOutput.TSestModel <- function(x)startOutput(x$data)

end.TSestModel <- function(x)end(x$data)
endInput.TSestModel <- function(x)endInput(x$data)
endOutput.TSestModel <- function(x)endOutput(x$data)

frequency.TSestModel <- function(x)frequency(x$data)
frequencyInput.TSestModel <- function(x)frequencyInput(x$data)
frequencyOutput.TSestModel <- function(x)frequencyOutput(x$data)

periods.TSestModel <- function(x)periods(x$data)
periodsInput.TSestModel <- function(x)periodsInput(x$data)
periodsOutput.TSestModel <- function(x)periodsOutput(x$data)

input.data.TSestModel <- function(x, series=seqN(nseriesInput(x)))
    input.data(x$data, series=series)
output.data.TSestModel <- function(x, series=seqN(nseriesOutput(x)))
    output.data(x$data, series=series)

nseriesInput.TSestModel <- function(x)  nseriesInput(x$data)
nseriesOutput.TSestModel <- function(x) nseriesOutput(x$data)

seriesNamesInput.TSestModel <- function(x)
 {m <- seriesNamesInput(x$model)
  d <- seriesNamesInput(x$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

seriesNamesOutput.TSestModel <- function(x)
 {m <- seriesNamesOutput(x$model)
  d <- seriesNamesOutput(x$data)
  if (!all(m==d)) 
    warning("data and model names do not correspond. Model names returned.")
  m
 }

seriesNames.TSestModel <- function(x)
  {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

"seriesNames<-.TSestModel" <- function(x, value)
   { seriesNamesInput(x) <- value$input
    seriesNamesOutput(x) <- value$output
    x
   }


"seriesNamesInput<-.TSestModel" <- function(x, value)
   {seriesNamesInput(x$model) <- value;
    seriesNamesInput(x$data ) <- value;
    x
   }
"seriesNamesOutput<-.TSestModel" <- function(x, value) 
   {seriesNamesOutput(x$model) <- value;
    seriesNamesOutput(x$data ) <- value;
    x
   }

identifiers.TSestModel <- function(obj){identifiers(TSdata(obj))}
sourcedb.TSestModel <- function(obj){sourcedb(TSdata(obj))}
sourceserver.TSestModel <- function(obj){sourceserver(TSdata(obj))}
source.info.TSestModel <- function(obj){source.info(TSdata(obj))}


############################################################################

#    methods for TSdata  <<<<<<<<<<
#      note there are some methods in dse1c too

############################################################################

is.TSdata <- function(obj) { inherits(obj, "TSdata")}

print.TSdata <- function(x, ...)
{  if(0 != (nseriesInput(x)))
     {cat("input data:\n")
      print(input.data(x),...)
      if(!is.null(x$input.transformations))
          {cat("input.transformations:\n")
           print(x$input.transformations, ...)
          }
      if(!is.null(seriesNamesInput(x)))
          {cat("input.names:\n")
           print(seriesNamesInput(x), ...)
          }
      cat("\n")
     }
  if(0 != (nseriesOutput(x)))
     {cat("output data:\n")
      print(output.data(x),...)
      if(!is.null(x$output.transformations))
         {cat("output.transformations:\n")
          print(x$output.transformations, ...)
         }
      if(!is.null(seriesNamesOutput(x)))
         {cat("output.names:\n")
          print(seriesNamesOutput(x), ...)
         }
     }
   cat("\n")
   if(!is.null(x$retrieval.date))
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) 
     {cat("source:\n")
      print(x$source)
     }
  invisible(x)
}

summary.TSdata <- function(object)
  {d <- output.data(object)
   if (!is.tframed(d)) d <- as.ts(d)
   st <- start(d)
   en <- end(d)
   fr <- frequency(d)
   d <- cbind(d,input.data(object))

   classed(list(  # summary.TSdata constructor
      description=object$description,
      start=st,
      end=en,
      freq=fr,
      sampleT= nrow(output.data(object)),
      p=nseriesOutput(object),
      m=nseriesInput(object),
      ave=apply(d,2,mean),
      max=apply(d,2,max),
      min=apply(d,2,min),
      retrieval.date=object$retrieval.date, 
      source=object$source), "summary.TSdata")
  }

print.summary.TSdata <- function(x, digits=options()$digits)
{  if (!is.null(x$description)) cat(x$description,"\n")
   cat("start.: ", x$start, " end.: ",  x$end," Frequency: ", x$freq,"\n")
   cat("sample length = ",x$sampleT,"\n")
   cat("number of outputs=",x$p, "   number of inputs=",x$m, "\n")
   cat("average :\n")
   print(x$ave)
   cat("max    :\n")
   print(x$max)
   cat("min    :\n")
   print(x$min)
   cat("\n")
   if(!is.null(x$retrieval.date)) 
      cat("retrieval date: ", x$retrieval.date, "   ")
   if(!is.null(x$source)) {cat("source:\n"); print(x$source)  }
   invisible(x)
}


tfplot.TSdata <- function(x, ..., start.=NULL,end.=NULL,
        select.inputs  = seq(length=nseriesInput(x)),
        select.outputs = seq(length=nseriesOutput(x)),
        Title=NULL, xlab=NULL, ylab=NULL, 
	graphs.per.page=5, mar=par()$mar, reset.screen=TRUE)
{# plot input and output data.
 # ... is a list of objects of class TSdata (with similar input
 #       and output dimensions.
 #previously mar = if(is.R()) c(3.1,6.1,1.1,2.1) else c(5.1,6.1,4.1,2.1)
 # Note that using ... like this means it cannot be used to pass additional
 #   arguments to plot, so unfortunately all necessary plot arguments must be 
 #   explicit in the arguments to tfplot.TSdata
 # start. is the starting point (date)  and end. the ending point for
 # plotting. If not specified the whole sample is plotted.
# output graphics can be paused between pages by setting par(ask=T).
  x <- freeze(x)
  Ngraphs <- length(select.outputs) + length(select.inputs)
  if(reset.screen)
    {Ngraphs <- min(Ngraphs, graphs.per.page)
     old.par <- par(mfcol = c(Ngraphs, 1), mar=mar) # previously c(5.1,6.1,4.1,2.1)) 
     on.exit(par(old.par))
    }
  if (is.null(Title)) Title <- ""
  if (!is.null(ylab))
       names <- rep(ylab, nseriesInput(x) + nseriesOutput(x))
  else names <-  c(seriesNamesInput(x),  seriesNamesOutput(x))
#  if (0 != length(select.inputs)) 
    {for (i in select.inputs) 
      {j <- 0
       z <- matrix(NA, periodsInput(x), length(append(list(x),list(...))))
       if(mode(i)=="character") i <- match(i, seriesNamesInput(x))
       for (d in append(list(x),list(...)) ) 
         {if (!is.TSdata(d))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
          d <- freeze(d)
          j <- j + 1
          z[, j] <- input.data(d, series = i)
         }
       tframe(z) <-tframe(input.data(x))
       if (!is.null(start.)) z <- tfwindow(z,start=start.)
       if (!is.null(end.))   z <- tfwindow(z,end=end.)
       tfOnePlot(z,ylab=names[i]) 
       if(i==1) title(main = Title)
    } }
  for (i in select.outputs) 
    {j <-0
     if(mode(i)=="character") i <- match(i, seriesNamesOutput(x))
     z <- matrix(NA,periodsOutput(x),length(append(list(x),list(...))))
     #z <-NULL 
     for (d in append(list(x),list(...)) ) 
       {if (!is.TSdata(d))
            stop("Expecting TSdata objects. Do not truncate argument names as that can cause a problem here.")
        d <- freeze(d)
        j <- j+1
        z[,j]<-output.data(d, series=i) 
       }
     tframe(z) <-tframe(output.data(x))
     if (!is.null(start.)) z <- tfwindow(z,start=start.)
     if (!is.null(end.))   z <- tfwindow(z,end=end.)
     tfOnePlot(z,ylab=names[nseriesInput(x) + i]) # tsplot
     if((0 == nseriesInput(x)) & (i==1)) title(main = Title)
    }
  invisible()
}
 
input.data.TSdata <- function(x, series=seqN(nseriesInput(x)))
  {if (is.null(x$input)) NULL  else select.series(x$input, series=series) }

output.data.TSdata <- function(x, series=seqN(nseriesOutput(x)))
  {if (is.null(x$output)) NULL  else select.series(x$output, series=series) }


"input.data<-.TSdata" <- function(x, value) 
   {cls <- dseclass(x); x$input <- value; dseclass(x) <- cls
    x
   }


"output.data<-.TSdata" <- function(x, value)
   {cls <- dseclass(x); x$output <-value; dseclass(x) <- cls
    x
   }

# Note: series names changed to an attribute of input and output for 
#     data but not for models!!!!

seriesNames.TSdata <- function(x)
 {list(input=seriesNamesInput(x), output=seriesNamesOutput(x))}

 seriesNamesInput.TSdata <- function(x) {seriesNames( input.data(x))}
seriesNamesOutput.TSdata <- function(x) {seriesNames(output.data(x))}

"seriesNames<-.TSdata" <- function(x, value)
   { seriesNamesInput(x) <-  value$input
    seriesNamesOutput(x) <-  value$output
    x
   }

"seriesNamesInput<-.TSdata" <- function(x,  value)
   {if (length( value) != nseriesInput(x))
       stop("number of series and number of names do not match.")
    attr(input.data(x), "seriesNames")  <- value
    x
   }
"seriesNamesOutput<-.TSdata" <- function(x,  value) 
    {if (length( value) != nseriesOutput(x))
       stop("number of series and number of names do not match.")
    attr(output.data(x), "seriesNames")  <- value
    x
   }

nseriesInput.TSdata <- function(x)
   {if (is.null(x$input)) 0 else nseries(x$input)}

nseriesOutput.TSdata <- function(x)
   {if (is.null(x$output)) 0 else nseries(x$output)}

start.TSdata <- function(x)
{i  <-  startInput(x)
 o  <- startOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 else return(c(i,o))
}

startInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(start(x$input))
 }

startOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(start(x$output))
 }

end.TSdata <- function(x)
{i  <-  endInput(x)
 o  <- endOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

endInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(end(x$input))
 }

endOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(end(x$output))
 }

frequency.TSdata <- function(x)
{i  <-  frequencyInput(x)
 o  <- frequencyOutput(x)
 if (((!is.null(o)) & (!is.null(i))) && all(i==o)) return(o)
 return(c(i,o))
}

frequencyInput.TSdata <- function(x)
 {if (is.null(x$input))  return(NULL)
  else return(frequency(x$input))
 }

frequencyOutput.TSdata <- function(x)
 {if (is.null(x$output))  return(NULL)
  else return(frequency(x$output))
 }


periods.TSdata <- function(x) UseMethod("periodsOutput")
periodsOutput.TSdata <- function(x)  dim(output.data(x))[1]
periodsInput.TSdata <- function(x)  dim(input.data(x))[1]

tbind.TSdata <- function(d1, d2)
 {if( ! (is.TSdata(d1) & is.TSdata(d2)))
     stop("tbind requires arguments to be of a similar type (ie. TSdata).")
  if ((0 != nseriesInput(d1)) || (0 != nseriesInput(d2)) )
    input.data(d1) <- tbind(input.data(d1),input.data(d2))
  if ((0 != nseriesOutput(d1)) || (0 != nseriesOutput(d2)) )
    output.data(d1) <- tbind(output.data(d1),output.data(d2))
  d1
 }



test.equal.TSdata <- function(obj1, obj2, fuzz=1e-16)
  {r <- TRUE
   if (r & (!is.null(obj1$input)))
     {if(is.null(obj2$input)) r <- FALSE
      else  r <-test.equal.matrix(obj1$input, obj2$input, fuzz=fuzz)
     }
   if (r & (!is.null(obj1$output)))
     {if(is.null(obj2$output)) r <- FALSE
      else r <-test.equal.matrix(obj1$output, obj2$output, fuzz=fuzz)
     }
   r
  }


identifiers.TSdata <- function(obj) {identifiers(obj$source)}
sourcedb.TSdata <- function(obj) {sourcedb(obj$source)}
sourceserver.TSdata <- function(obj) {sourceserver(obj$source)}
source.info.TSdata <- function(obj) {source.info(obj$source)}




TSdata <- function (data=NULL, ...) {UseMethod("TSdata")} 
 

TSdata.default <- function(data=NULL, input=NULL, output=NULL, ...)  
{if (is.null(data) && (!is.null(input) | !is.null(output) ))
    {if(!is.null(input) && is.vector(input))
           input <- matrix(input, length(input),1)
     if(!is.null(output) && is.vector(output))
           input <- matrix(output, length(output),1)
     data <- classed(list(input=input, output=output), "TSdata") # constructor
  }else 
     data <- classed(data, "TSdata")   # constructor keeps other list elements
  
 if(!is.list(data)) stop("TSdata.default could not construct a TSdata format.")	    
 if   ( 0 == nseriesInput(data)   &    0 == nseriesOutput(data) 
    |  (0 != nseriesOutput(data) && !is.matrix(output.data(data)))
    |  (0 !=  nseriesInput(data) && !is.matrix( input.data(data))) )
      stop("TSdata.default could not construct a TSdata format.")
 data
}

TSdata.TSdata <- function(data, ...)  {data} # already TSdata 

TSdata.TSestModel <- function(data, ...)  {data$data} # extract TSdata 

as.TSdata <- function(d) 
 {# Use whatever is actually $input and $output,
  #strip any other class, other parts of the list, and DO NOT use
  # input.data(d) and output.data(d) which may eg reconstitute something
  TSdata(input=d$input, output=d$output)
 }

"tframe<-.TSdata" <- function(x, value) {
 if (0 != nseriesInput(x))  tframe(input.data(x))  <- value
 if (0 != nseriesOutput(x)) tframe(output.data(x)) <- value
 x
}

tframed.TSdata <- function(x, tf=NULL, input.names=NULL, output.names=NULL)  
{# switch to tframe representation
 if(0 != (nseriesOutput(x)))
       output.data(x) <-tframed(output.data(x), tf=tf, names=output.names)
 if (0 != (nseriesInput(x)))
        input.data(x) <-tframed(input.data(x),  tf=tf, names= input.names)
 x
}  


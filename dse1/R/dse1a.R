
# For installation instructions see the file read.me or the brief user's
#    guide (postscipt file guide.ps).

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


   .DSECOMPILED <- T

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
      ok <-      require("syskern", warn.conflicts=F)
      ok <- ok & require("tframe",  warn.conflicts=F)
      if(!ok) warning("This package requires the syskern and tframe packages.")


      global.assign("DSE.HOME", DSE.HOME)            
#      if (pmatch("package:dse1", search()) >4)
#         warning("The DSE library may not work properly if it is not near the beginning of the search list.")
#     source(paste(DSE.HOME,"/data/egJofF.1dec93.data.R", sep=""))
#     eg1.DSE.data <- 
#           example.get.eg.raw.data(paste(DSE.HOME,"/data/eg1.dat", sep=""))

# in very old Guide eg1.DSE.data       was called  eg1.DSE.data.diff.all.raw
# previously        eg1.DSE.data.diff  was called  eg1.DSE.data.diff.all

#     take the log and difference of the output variables and puts the  
#     data in the variable eg1.DSE.data.diff.
#  global.assign("eg1.DSE.data.diff", example.convert.eg.data(eg1.DSE.data) )

#  example.BOC.93.4.data.trunc <- example.truncate.data(eg1.DSE.data.diff)

     } else
   if (is.S())
     {global.assign("DSE.HOME", DSE.HOME)
      if (pmatch(DSE.HOME,search()) >2)
         warning("The DSE library may not work properly if it is not near the beginning of the search list. Use library(..., first=T)")
     }

   #paste("Setting DSE.HOME to ", DSE.HOME)

# Once these functions and data have been loaded then an example 
# can be run by:
#        example.VAR.SVD(example.BOC.93.4.data.trunc)
# also recommended is:
# example.verify.data(eg1.DSE.data.diff) # print some summary statistics

   invisible( if (.DSECOMPILED) load.DSE.fortran() else T )
  }



# The following variable is used in Windows dll.load
dse.win.for.tab <- c("simss" ,"smooth" ,"kfp" ,"kfprj" ,"kfepr" ,"kf" ,"simarma"
          ,"armap" ,"armaprj" ,"armaepr" ,"arma" ,"dataepr","inverse"
          ,"gend","cstat_f","efcurve_f","rlcurve_f","wep")

load.DSE.fortran <- function(from=paste(DSE.HOME,"/lib", sep=""), large=F)
{#   if(is.MSwindows()) from <- paste(from, "/_Data", sep="")
 #   else               from <- paste(from, "/.Data", sep="")
 if ( (0 != nchar(from)) &
     "/" != substring(from, first=nchar(from))) from <- paste(from,"/", sep="")
 if (is.R()) r <- library.dynam("dse1")  # does not use from
 else if (is.S())
    {if(is.MSwindows()) 
          {r <- dll.load(paste(from,"dsefor.dll", sep=""),dse.win.for.tab)
           if (1!=r) warning("dll load was NOT successful.")
          }
     else 
        {if (large) r <- dyn.load(paste(from,"dsefor.large.o", sep=""))
         else       r <- dyn.load(paste(from,"dsefor.o", sep=""))
    }   }
 invisible(r)
} 



help.start.DSE <- function (browser = "netscape") 
   {#browser = "netscape -mono"
    #browser = "mosaic"
    #redirecting standard input and output is necessary in order to return to S 
    system.call(paste(browser, " ", DSE.HOME,"/help/dsehome.htm </dev/null >/dev/null &", sep=""))
    invisible()
   }

##############################################################################

#  end of section containing installation dependent strings.

##############################################################################

##############################################################################

#  section containing documentation "stubs" (specific methods 
#    for generic functions) and utility functions
#       so that R CMD build does not complain.

##############################################################################










       


##############################################################################

#  end of section containing documentation "stubs" (specific methods 
#    for generic functions) and utility functions
#       so that R CMD build does not complain.

##############################################################################

# Some utilities for insulating DSE code from S/R changes.

# See also classed, tfclass, and tfclassed in tframe.

##############################################################################

if (is.R())
  { dseclass    <- .Alias(class)
   "dseclass<-" <- .Alias(get("class<-"))
    dsescan <- function(file="",quiet=T, ...){scan(file=file, quiet=quiet, ...)}
  } else
if (is.S())
  { dseclass    <- class
   "dseclass<-" <- function(x, value){ class(x) <- value ; x }
    dsescan <- function(file="",quiet=T, ...){scan(file=file, ...)}
  }


##############################################################################

##############################################################################



# The following two functions were used experimentally but have been disabled
#    because of the possibility they will cause problems in functions
#    outside this library.


# "%$%" <- function(x, n)
#   {#  This function allows %$% to be used in place of $ to require exact
#    #     matching of the name in the list.
#    #  However, beware that the order of precedence is not the same as $ and
#    #    may thus cause problems.
#    #  Alternately use "$<-" <- get("%$%<-") and "$" <- get("%$%")
#    #    to replace default action of $,  but beware.
#    if ( is.null(x)) return(NULL)
#    if (!is.list(x)) stop("LHS of $ is not a list")
#    n <- as.character(sys.call())[3]
#    if ( is.null(n)) return(NULL)
#    n <- match(n, names(x))
#    if ( is.na(n)) return(NULL)
#    x[[n]]
#   }

# "%$%<-" <- function(x,n, value)
#   {if ( is.null(x)) x <-list()
#    if (!is.list(x)) stop("assignment attempted to object which is not a list")
#    n <- as.character(sys.call())[3]
#    x[[n]] <- value 
#   x
#   }

      # the following would enforce exact argument matching in lists, which 
      #  seems preferable, but it is not used because it could potentially cause
      #  problems if users take  advantage of truncated matching in their 
      #  functions or use.
      #"$" <- get("%$%")  # enforce exact matching !!! do in .First.lib and warn
      #"$<-" <- get("%$%<-")
      

##############################################################################



version.dse <- function() cat(dse.version.information, "\n")

##############################################################################

# The code in files dse1a.s, dse1b.s dse1c.s and dse1d.s was divided roughly into 
#   the groups listed below, but the organization has changed a little bit.
#   The actual grouping can be seen by grep ing on the string '<<<<<<' eg:
#    grep "<<<<<<" dse1*.s

#    Functions which work on a model (i.e. if a model with data is allowed as
#             an arguement then the data is ignored):
#        -model summary, description, display, comparison and 
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
#        -data interface example functions

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

#     functions for model summary, description, display and
#      comparison and functions for calculation of model properties

############################################################

display <- function(x,...)  UseMethod("print")

print.TSestModel <- function(x, ...) 
{ cat("neg. log likelihood=",x$estimates$like[1],"\n")
  if(!is.null(x$converged)) {if(!x$converged) cat(" (not converged)") }
  print(x$model,...) 
  invisible(x)
}

print.SS <- function(x, digits=4) 
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

print.ARMA <- function(x, digits=4, L=T, fuzz=1e-10) 
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
   residual <- object$estimates$pred[,,drop=F] - output.data(object)[,,drop=F]
   sampleT <- nrow(residual)
   p <- ncol(residual)	
   Om <- t(residual) %*% residual/sampleT
   rmse <- matrix( diag(Om)^.5 ,1,p)
   dimnames(rmse) <- list(c("RMSE"), output.series.names(object))

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
  {m <- input.dimension(object)
   p <- output.dimension(object)
   n <- nrow(object$F)
   classed(list(  # summary.SS constructor
         description=object$description,
         input.series=input.series.names(object),
         output.series=output.series.names(object),
         innov=is.innov.SS(object),
         m=m,
         p=p,
         n=n,
         P=n * (m+2*p),  #assumes full rank noise
         P.actual=length(object$parms),
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

summary.ARMA <- function(object,digits=6)
  {m <- input.dimension(object)
   p <- output.dimension(object)
   classed(list(  # summary.ARMA constructor
         description=object$description,
         input.series=input.series.names(object),
         output.series=output.series.names(object),

         a=dim(object$A)[1]-1,
         b=dim(object$B)[1]-1,
         c=dim(object$C)[1]-1, 
         m=m,
         p=p,
         P.actual=length(object$parms),
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
 


tfplot.TSmodel <- function(x, ...)
 {stop("This is just a TSmodel. plot needs a TSestModel or data.") }


tfplot.TSestModel <- function(..., start.=NULL,end.=NULL,Title=NULL, 
  reset.screen=T, select.inputs=NULL, select.outputs=NULL, graphs.per.page=5)
{
# plot one-step ahead estimates and actual data.
# ... is a list of models of class TSestModel.
# start. is the starting point (date) for plotting.
  model <-list(...)[[1]]
  if (!is.TSestModel(model)) TS.error.exit(clss="TSestModel") 
  if (is.null(Title))
     Title <- "One step ahead predictions (dotted) and actual data (solid)"
  p<- output.dimension(model)
  if (is.null(select.outputs)) select.outputs <-1:p
  if (all(0==select.outputs)) select.outputs <- NULL
  Ngraphs <- length(select.outputs)
  if (!is.null(select.inputs)) if (all(0==select.inputs)) select.inputs <- NULL
  m<- input.dimension(model)
  if (is.null(m)) m <-0
  else Ngraphs <- Ngraphs+length(select.inputs)  # NULL is zero
  Ngraphs <- min(Ngraphs, graphs.per.page)
  if(reset.screen) 
    {old.par <- par(mfcol = c(Ngraphs, 1), mar= c(5.1,6.1,4.1,2.1))
     on.exit(par(old.par))
    }
  names <-output.series.names(model)
  if (m!=0) names <-c(input.series.names(model), names)
  if (is.null(names)) names <- rep(" ",m+p)
  if (m!=0) 
    {for (i in select.inputs) 
      {z <-NULL 
       for (model in list(...) ) z<-tbind(z,input.data(model, series=i))
       tframe(z) <-tframe(input.data(model))
       if (is.null(start.)) start.<-start(z)
       if (is.null(end.))   end.  <-end(z)
       tfplot(tfwindow(z,start=start.,end=end., warn=F),ylab=names[i]) # tsplot
       if(i==1) title(main = Title)
    } }
  for (i in select.outputs ) 
    {z <-c(output.data(model, series=i),
           rep(NA,dim(model$estimates$pred)[1]-periods(model$data)))
     for (model in list(...)) z<-cbind(z,model$estimates$pred[,i,drop=F])
     tframe(z) <- tframe(output.data(model))
     if (is.null(start.)) start.<-start(z)
     if (is.null(end.))   end.  <-end(z)
     tfplot(tfwindow(z,start=start.,end=end., warn=F),ylab=names[m+i]) # tsplot
     if(i==1) title(main = Title)
    }
  invisible()
}
    


test.equal.TSestModel <- function(obj1, obj2, ...) # this could be better
  { test.equal.TSmodel( obj1$model, obj2$model, ...) &
        test.equal.TSdata(obj1$data, obj2$data, ...)
  }

test.equal.TSmodel<- function (model1,model2, fuzz=0)
{# return T if models are identical (excluding description)
  r       <- all(dseclass(model1) == dseclass(model2))
  if (r) r <-length(model1$parms) == length(model2$parms)
  if (r) r <-all(fuzz >= abs(model1$parms   -     model2$parms))
  if (r) r <-length(model1$location) == length(model2$location)
  if (r) r <-all(model1$location  ==     model2$location)
  if (r) r <-length(model1$i) == length(model2$i)
  if (r) r <-all(model1$i  ==     model2$i)
  if (r) r <-length(model1$j) == length(model2$j)
  if (r) r <-all(model1$j  ==     model2$j)
  if (r) r <-length(model1$const) == length(model2$const)
  if (r) r <-all(model1$const  ==     model2$const)
  if (r) r <-length(model1$const.location) == length(model2$const.location)
  if (r) r <-all(model1$const.location  ==     model2$const.location)
  if (r) r <-length(model1$const.i) == length(model2$const.i)
  if (r) r <-all(model1$const.i  ==     model2$const.i)
  if (r) r <-length(model1$const.j) == length(model2$const.j)
  if (r) r <-all(model1$const.j  ==     model2$const.j)
  if (r)
    {if (is.ARMA(model1))    r <-test.equal.ARMA(model1,model2, fuzz=fuzz)
     else if (is.SS(model1)) r <-test.equal.SS(model1,model2, fuzz=fuzz)
    }
  r
}

test.equal.ARMA<- function (model1,model2, fuzz=0)
{    r <-length(model1$l) == length(model2$l)
     if (r) r <-all(model1$l  ==     model2$l)
     if (r) r <-length(model1$const.l) == length(model2$const.l)
     if (r) r <-all(model1$const.l  == model2$const.l)
     if (r) r <-length(model1$A) == length(model2$A)
     if (r) r <-all(fuzz >= abs(model1$A   -     model2$A))
     if (r) r <-length(model1$B) == length(model2$B)
     if (r) r <-all(fuzz >= abs(model1$B   -     model2$B))
     if (r) 
           {if (is.null(model1$C)) r <- is.null(model2$C)
             else
               {if (r) r <-length(model1$C) == length(model2$C)
                if (r) r <-all(fuzz >= abs(model1$C   -     model2$C))
           }  }
  r
}

test.equal.SS<- function (model1,model2, fuzz=0)
{    r <-length(model1$F) == length(model2$F)
     if (r) r <-all(fuzz >= abs(model1$F   -     model2$F))
     if (r) 
         {if(is.null(model1$G)) r <- is.null(model2$G)
          else
            {if (r) r <-length(model1$G) == length(model2$G)
             if (r) r <-all(fuzz >= abs(model1$G   -     model2$G))
         }  }
     if (r) r <-length(model1$H) == length(model2$H)
     if (r) r <-all(fuzz >= abs(model1$H   -     model2$H))
     if (is.innov.SS(model1))
       {if (r) r <-length(model1$K) == length(model2$K)
        if (r) r <-all(fuzz >= abs(model1$K   -     model2$K))
       }
     else
       {if (r) r <-length(model1$Q) == length(model2$Q)
        if (r & (0 != length(model2$Q)) ) r <-all(fuzz >= abs(model1$Q - model2$Q))
        if (r) r <-length(model1$R) == length(model2$R)
        if (r & (0 != length(model2$R)) ) r <-all(fuzz >= abs(model1$R - model2$R))
       }
     if (r) if(is.null(model1$z0)) r <- is.null(model2$z0)
     else
       {if (r) r <-length(model1$z0) == length(model2$z0)
        if (r) r <-all(fuzz >= abs(model1$z0   -     model2$z0))
       }
     if (r) if(is.null(model1$P0)) r <- is.null(model2$P0)
     else
       {if (r) r <-length(model1$P0) == length(model2$P0)
        if (r) r <-all(fuzz >= abs(model1$P0   -     model2$P0))
       }
  r
}


McMillan.degree <- function (model,  ...) UseMethod("McMillan.degree")

McMillan.degree.TSestModel <- function(model, ...)
 {McMillan.degree(TSmodel(model),...) }

McMillan.degree.ARMA<- function   (model, fuzz=1e-4, verbose=T, warn=T) {
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

McMillan.degree.SS<- function  (model, fuzz=1e-4) {
   cat("state dimension = ",dim(model$F)[1],"/n")
   invisible()
}

stability <- function(obj, ...) UseMethod("stability")

stability.TSestModel <- function(model, ...){stability(TSmodel(model),...)}
stability.TSmodel <- function(model, fuzz=1e-4, ...)
  {stability(roots(model, fuzz=fuzz, randomize=F), ...)}

stability.roots  <- function (z, digits=8, verbose=T) 
   {#z <- roots(model, fuzz=fuzz, randomize=F)
    if (all(Mod(z) < 1.0)) s <- T
    else                   s <- F
    if (verbose)
      {cat("Eigenvalues of F and moduli are:\n")
       print(cbind(z,Mod(z)),digits=digits)
       if (s) cat("The system is stable.\n")
       else   cat("The system is NOT stable.\n")
      }
    s
   }

stability.ARMA  <- function (model, fuzz=1e-4, digits=8, verbose=T) 
   {z <- roots(model, fuzz=fuzz, randomize=F)
    if (all(Mod(z) < 1.0)) s <- T
    else                   s <- F
    if (verbose)
      {cat("Distinct roots of det(A(L)) and moduli are:\n")
       print(cbind(1/z,Mod(1/z)),digits=digits)
       cat("\nInverse of distinct roots of det(A(L)) and moduli are:\n")
       print(cbind(z,Mod(z)),digits=digits)
       if(!is.null(model$TREND)) cat("Trend not taken into account: ")
       if (s) cat("The system is stable.\n")
       else   cat("The system is NOT stable.\n")
      }
    s
   }


roots <- function(obj, ...) UseMethod("roots")

roots.TSestModel <- function(model, ...){roots(TSmodel(model),...)}

roots.SS  <- function (model, fuzz=0, randomize=F) 
{   z <- eigen(model$F)$values
    if (randomize) if (sample(c(T,F),1)) z <- Conj(z)
      #this prevents + - ordering of complex roots (for Monte Carlo evaluations)
    classed(z,"roots")  # constructor (roots.SS)
}


roots.ARMA  <- function (model, fuzz=0, verbose=T, randomize=F, warn=T, by.poly=F) 
{   if(by.poly) z <- 1/polyroot.det(model$A)
    else        z <- roots(to.SS(model))
    if (fuzz!=0)
      {zz <- outer(z,z,FUN="-")  # find distinct roots within fuzz
       z <- z[ !apply((outer(1:length(z),1:length(z),FUN="<")
                 & (Mod(zz) <fuzz)),2,any)]
      }
    # add unit roots for TREND elements.
    if (!is.null(model$TREND))
      {z <- c(rep(1,sum(0!=model$TREND)), z)
       if (warn)
         warning("Unit roots have been added for non-zero trend elements.")
      }
    if (randomize) if (sample(c(T,F),1)) z <- Conj(z)
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
observability.TSestModel <- function(Emodel){observability(TSmodel(model)) }

observability.SS  <- function (model)
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
observability.ARMA  <- function (model)
{ cat("not applicable to ARMA models\n")
  invisible()
}


reachability<- function (model)  
 {# calculate singular values of reachability matrix
  UseMethod("reachability")
}

reachability.TSestModel <- function(Emodel){reachability(TSmodel(model))}

reachability.SS<- function (model)
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

reachability.ARMA<- function (model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


check.balance <- function(model)  
 { # calculate the difference between observability and reachability gramians 
  UseMethod("check.balance")
 }
check.balance.TSestModel <- function(model){check.balance(TSmodel(model))}

check.balance.SS  <- function (model){ 
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
check.balance.ARMA  <- function  (model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


check.balance.Mittnik<- function (model)  
{# calculate the difference between observability and controllability
 #   gramians with model transformed a la Mittnik 
 UseMethod("check.balance.Mittnik")
}

check.balance.Mittnik.TSestModel <- function(model)
   {check.balance.Mittnik(TSmodel(model))}

check.balance.Mittnik.SS  <- function (model){ 
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
check.balance.Mittnik.ARMA  <- function  (model){ 
  cat("not applicable to ARMA models\n")
  invisible()
}


############################################################

#     functions for model conversion   <<<<<<<<<<

############################################################

to.SS <- function(model, ...) {UseMethod("to.SS") }
to.SS.TSestModel <- function(model, ...) 
	{l(to.SS(TSmodel(model), ...),TSdata(model))}
to.SS.SS <- function(model, ...) {model}

to.SS.ARMA <- function(model,...)
{# convert an ARMA (or VAR) to a SS (innovations) representation
    if (is.null(model$A)) a<-0
    else a <- dim(model$A)[1] - 1  #order of polynomial arrays
    if (is.null(model$B)) b<-0
    else b <- dim(model$B)[1] - 1
    if (is.null(model$C)) cc<-0
    else cc<- dim(model$C)[1] - 1
    if ((b<=a) & (cc<=(a-1))) model <- to.SS.augment(model)
    else  model <-to.SS.nested(model,...) #  (otherwise best working method) 
                  # A better approach would be an algorithm like Guidorzi's. 
 model
}


to.SS.nested <- function(model, ...) {UseMethod("to.SS.nested") }
to.SS.nested.TSestModel <- function(model, ...)
   {to.SS.nested(TSmodel(model), ...)}

to.SS.nested.SS <- function (model, n=NULL, Aoki=F)
{# convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <-ncol(model$F)  
  if (Aoki) return(Aoki.balance(model, n=n))
  else      return(balance.Mittnik(model, n=n)) 
}

to.SS.nested.ARMA  <- function (model, n=NULL, Aoki=F)
{# convert to a nested-balanced state space model by svd  a la Mittnik (or Aoki)
  if (is.null(n)) n <- McMillan.degree.calculation(model)$distinct
  if (Aoki) return(Aoki.balance(model, n=n))
  else      return(balance.Mittnik(model, n=n)) 
}


to.SS.augment <- function(model, ...) {UseMethod("to.SS.augment") }
to.SS.augment.TSestModel <- function(model, ...)
   {l(to.SS.augment(TSmodel(model), ...), TSdata(model))}


to.SS.augment.ARMA <- function(model, fuzz=1e-14) 
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
   p <- output.dimension(model)          # Dim of endoenous Variables.
   m <-  input.dimension(model)          # Dim of exogenous Variables.
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
         input.names= input.series.names(model),
        output.names=output.series.names(model))        
 }


gmap <- function(g, model) 
{# convert to an equivalent representation using a given matrix
 if(!is.TSm.or.em(model)) TS.error.exit()
 if (is.TSestModel(model)) model <- model$model
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
  if(!is.TSm.or.em(model1)) TS.error.exit()
  if (is.TSestModel(model1)[1]) model1 <- model1$model
  if(!is.TSm.or.em(model2)) TS.error.exit()
  if (is.TSestModel(model2)[1]) model2 <- model2$model
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
 
  if(!is.TSm.or.em(model)) TS.error.exit()
  if (is.TSestModel(model)) model <- model$model
  if (is.null(constants))
    {p <-abs(model$parms-1.0) < fuzz
     model$const <- c(model$const,rep(1.0,sum(p)))
     model$const.location <- c(model$const.location,model$location[p])
     model$const.i <- c(model$const.i,model$i[p])
     model$const.j <- c(model$const.j,model$j[p])
     if(is.ARMA(model)) model$const.l <- c(model$const.l,model$l[p])
     p <- (!p) & (abs(model$parms) > fuzz) 
     model$parms <-model$parms[p]
     model$location <- model$location[p]
     model$i <- model$i[p]
     model$j <- model$j[p]
     if(is.ARMA(model)) model$l <- model$l[p]
     return(set.arrays(model))
    }
  else
    return(set.parameters(TSmodel(append(model,constants))))
}


to.parsim <- function(model, ...) {UseMethod("to.parsim") }
to.parsim.TSestModel <- function(model, ...)
   {l(to.parsim(TSmodel(model), ...), TSdata(model))}


to.parsim.TSmodel <- function(model, minf = nlmin, small=1e-5, equiv=F,
                   max.iter=100,max.fcal=200)
{ # find a model with fewer non-constant matrix entries.. not working well
# Find a parsimonious almost equivalent ARMA model by finding an invertible
# pxp matrix which pre-multiplies A,B, and C to make more entries zero or one
# (1 & 0 are treated as constants not as parameters).
# Numerically it is necessary to set very small parameters to 0 (.. 1.0). 
# The judgement about small is crude.
# The idea in the objective func is to highly penalize small deviations from
# 1 or zero, but not penalize large deviations much more.
# If equiv is TRUE then the resulting model is forced to have the same noise cov
# determinant (up to numerical accuracy, which can be bad if the system is
#  degenerate and det(cov) is very small), otherwise,  parameters may be further 
#  eliminated but the cov of the residual and
# the resulting likelihood (|cov| component) will change.
# This may actually be a linear problem? Using least squares would be faster but 
#  requires some more work!
  warning("This function does not work yet.")
  if ( !is.ARMA(model)) stop("parsim still only attempted for ARMA models")
	A <- model$A
 	B <- model$B
	C <- model$C
	global.assign("A",A,frame=1)
	global.assign("B",B,frame=1)
	global.assign("C",C,frame=1)
	global.assign("equiv",equiv,frame=1)
	funcAC <- function(para)
	{
		zA <- A	
		zB <- B 
		zC <- C
		p <- dim(zA)[2]
		g <- matrix(para, p, p)
		g[p, p] <- p - sum(diag(g)[1:(p - 1)])	# force trace =p
		for(l in 1:dim(zA)[1]) zA[l,  ,  ] <- g %*% A[l,  ,  ]	
#		for(l in 1:dim(zB)[1]) zB[l,  ,  ] <- g %*% B[l,  ,  ]
		for(l in 1:dim(zC)[1]) zC[l,  ,  ] <- g %*% C[l,  ,  ]
		error <- c(zA, zC)
#		sum(log(1e-100+pmin(error^2, (error-1)^2))) # may not be differentiable.
#		sum(log(1e-100+error^2))
		sum(error^2)
	}
	funcB <- function(para)
	{
		zB <-B # opt separately by post mult.(affects only cov of Om)
		p <- dim(zB)[2]
		g <- matrix(para, p, p)
		g[p, p] <- p - sum(diag(g)[1:(p - 1)])	 # force trace =p
		if (equiv) g <-g/(prod(svd(g)$d)^(1/p)) # force det =1
		for(l in 1:dim(zB)[1]) zB[l,  ,  ] <- B[l,  ,  ] %*% g
		error <- c(zB)
#		sum(log(1e-100+pmin(error^2, (error-1)^2)))
#		sum(log(1e-100+error^2))
		sum(error^2)
	}
	p <- dim(A)[2]               # first A and C
	para <- c(diag(1, p))
	para <- para[1:(p^2 - 1)]
	para <- minf(funcAC, para, max.iter=max.iter, max.fcal=max.fcal)
	print(para[[2]])
	print(para[[3]])
	g <- matrix(para[[1]], p, p)
	g[p, p] <- p - sum(diag(g)[1:(p - 1)])
	for(l in 1:dim(A)[1]) A[l,  ,  ] <- g %*% A[l, ,]	
	for(l in 1:dim(B)[1]) B[l,  ,  ] <- g %*% B[l, ,]
	for(l in 1:dim(C)[1]) C[l,  ,  ] <- g %*% C[l, ,]
	A[abs(A) < small] <- 0.0
	B[abs(B) < small] <- 0.0
	C[abs(C) < small] <- 0.0
	A[abs(A-1) < small] <- 1.0
	B[abs(B-1) < small] <- 1.0
	C[abs(C-1) < small] <- 1.0
	para <- c(diag(1, p))           # now B
        para <- para[1:(p^2 - 1)]
	para <- minf(funcB, para, max.iter=max.iter, max.fcal=max.fcal)
	print(para[[2]])
	print(para[[3]])
	g <- matrix(para[[1]], p, p)
	g[p, p] <- p - sum(diag(g)[1:(p - 1)])
        if (equiv) g <-g/prod(svd(g)$d) # force det =1
	for(l in 1:dim(B)[1]) B[l,  ,  ] <- B[l, ,] %*% g
	B[abs(B) < small] <- 0.0
        B[abs(B-1) < small] <- 1.0
        model$A <- A
	model$B <- B
	model$C <- C
	set.parameters(model)
}


to.SS.innov <- function(model)
{ # convert to an equivalent state space innovations representation
# This assumes that the noise processes in the arbitrary SS representation are 
#   white and uncorrelated.
 if(!is.TSm.or.em(model)) TS.error.exit()
 if (is.TSestModel(model)) model <- model$model
 if (!is.SS(model))  model <- to.SS(model)
 if ( is.non.innov.SS(model)) 
   {model$K <- model$Q %*% solve(model$R)
    model$R <- NULL
    model$Q <- NULL
   }  
 classed(model, c("innov","SS","TSmodel"))  # bypass constructor
}



to.SS.Oform <- function(model, ...) {UseMethod("to.SS.Oform") }

to.SS.Oform.TSestModel <- function(model) 
  {l(to.SS.Oform(TSmodel(model)), TSdata(model))
  }

to.SS.Oform.TSmodel <- function(model)
{# convert to a SS innovations representation with a minimum number 
# of parameters by converting as much of H as possible to I matrix.
# Any remaining reductions are done by converting part of ?? to I.
# It seems there should remain n(m+2p) free parameters in F,G,H,K, and Om is 
#  determined implicitly by the residual.
 if (!is.SS(model))       model <- to.SS(model)
 if (!is.innov.SS(model)) model <- to.SS.innov(model)
 n <- dim(model$H)[2]
 p <- output.dimension(model)
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
    cat("This procedure is not working yet.") 
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
  if(!is.TSm.or.em(model)) TS.error.exit()
  if (is.TSestModel(model)) model <- model$model
  if (!is.SS(model))         model <- to.SS(model)
  if ( is.non.innov.SS(model))  model <- to.SS.innov(model)
  p <-model$location == "f"
  model$const <- c(model$const,model$parms[p])
  model$const.location <- c(model$const.location,model$location[p])
  model$const.i <- c(model$const.i,model$i[p])
  model$const.j <- c(model$const.j,model$j[p])
  p <- !p
  model$parms <-model$parms[p]
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


to.SS.Chol <- function(model, ...) {UseMethod("to.SS.Chol") }

to.SS.Chol.TSestModel <- function(model, Om=NULL) 
  {if(is.null(Om)) Om <-model$estimates$cov
   l(to.SS.Chol(TSmodel(model), Om=Om), TSdata(model))
  }

to.SS.Chol.TSmodel <- function(model, Om=diag(1,output.dimension(model)))
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


to.ARMA <- function(model) {UseMethod("to.ARMA") }

to.ARMA.TSestModel <- function(model, ...) 
	{l(to.ARMA(TSmodel(model), ...), TSdata(model))}

to.ARMA.ARMA <- function(model) {model}

to.ARMA.SS <- function(model)
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
 ARMA(A=A,B=B,C=C, input.names =  input.series.names(model),
                  output.names = output.series.names(model))
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

acf.M.TSdata <- function(data, lag=round(6*log(periods(data))), 
           type ="covariance", sub.mean=T)
{#Estimate auto covariances from data and return as first block row of Hankel.
 # i.e. a matrix with partitions [M0|...|Mi|...|Ml] giving the cov calculated from 
 #  the data, including the exogenous series and return as first block row of Hankel.
 #  Each Mi is a p by (m+p) matrix, (m is the dimension of the exogeneous 
 #  series and p is the dimension of endogeneous series)
 #  ie.  Mi = [ cov{y(t),y(t-i)} | cov{y(t),u(t-i)} ] 
 #  N.B. - The part of the first block corresponding to y(t)y(t) may need to be discarded
 #     for Aoki's technique.  This will reverse the order of y and u!!!!
 # if type == "correlation" the result is scaled to give autocorrelations.
  data <- freeze(data)
  p <- ncol(output.data(data))
  m <- input.dimension(data)
  if (is.null(m)) m <-0
  d <- cbind(output.data(data),input.data(data))
  sampleT <-periods(data)
  if (sub.mean) d <- d- t(matrix(apply(d,2,mean),dim(d)[2],sampleT))
#  if (type == "correlation") d <- d %*% diag(1/diag(t(d)%*%d/sampleT)^.5)
#  z <-acf(d, type = type, plot=F)$acf
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

acf.M.TSestModel <- function(model, ...) {acf.M(TSmodel(model), ...)}

acf.M.TSmodel <- function(model, lag=NULL, type ="covariance", Psi=NULL)
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
 
if (!is.SS(model)) model <- to.SS(model)
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



Riccati <- function(A, B, fuzz=1e-10, iterative=F)
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
{
if(!is.TSm.or.em(model)) TS.error.exit()
if (is.TSestModel(model)) model <- TSmodel(model)
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
    stop <- F
    while((i<=blocks) & !stop)   #  no. of blocks affect Hankel size
       {FnGK <- FF %*% FnGK
        M <- cbind(M,H %*% FnGK)
        i <-i+1   # count should not be necessary, but insures an end.
        stop <- (i>3) & ( max(abs(FnGK)) < 1e-15)
       }
  }
else stop("markov.parms requires an ARMA or SS model.")
M       
}



############################################################

#     polynomial utility functions   <<<<<<<<<<

############################################################


old.polyprod <- function(a,b)
{ # product of two polynomials.
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
   if (is.vector(b) && (is.array(a) | is.matrix(a)))
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
               for (i in 1:(dim(a)[1])) 
                  for (j in 1:(dim(b)[3]))
                    for (k in 1:(dim(a)[2]))
                       r[,i,j] <- psum(r[,i,j], pprod(a[i,k],b[,k,j]))
              }
          }
   else if (is.array(a))
        if (is.matrix(b))
          {if (dim(a)[3] != dim(b)[1])
              stop("Matrix polynomial dimensions do not conform.")
           r <- array(0,c(dim(a)[1],dim(a)[2], dim(b)[2]))
           for (i in 1:(dim(a)[2])) 
              for (j in 1:(dim(b)[2]))
                 for (k in 1:(dim(a)[3]))
                    r[,i,j] <- psum(r[,i,j], pprod(a[,i,k],b[k,j]))
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
                 #array above insures b a matrix, drop=F cannot be used
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
                 #array above insures b a matrix, drop=F cannot be used
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



polysum        <- function # sum of two polynomials (including polynomial arrays)
   (a,b){
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

polyroot.det<- function (A)
{#roots of the determinant of A.  Note: polydet is slow. There is room for improvement here!
z <- polydet(A)
if (length(z)==1) 
  stop(paste("root cannot be calculated for degree 0 determinant = ", as.character(z)))
polyroot(z)
}

polydet <- function (A)
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
          polyprod(A[,i,1],Recall(A[,(1:n)[i!=(1:n)],2:n,drop=F])))
              r[is.na(r)] <- 0
              if (any(0==r)) 
                   {if (all(r==0)) r <- 0
                    else  r <-r[0 == rev(cumprod(rev(r==0)))] #remove trailing zeros
                    if(0==length(r)) r <- 0
             }     } 
          }
r
}

poly.value     <- function (coef,z)
{# evaluate a polynomial given by coef (constant first) at z
 # could be extended for matrix coef.
#           n-1           n-2
#  coef[n]*z   + coef[n-1]*z   + ... + coef[1]  
  coef %*% z^(0:(length(coef)-1))
}


characteristic.poly<- function # coefficients of the characteristic polynomial of a matrix
   (A) {
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
   {err <-T
    while (err)
     {cat(prmt)
      n <- as.integer(readline()) # crude. this truncates reals
      if (is.na(n)) cat("value must be an integer\n")
      else err <- F
     }
    n
   }




input.dimension <- function(x, ...)UseMethod( "input.dimension")
input.dimension.SS <- function(obj)
   {if (is.null(obj$G)) return(0)
    else   return(dim(obj$G)[2])
   }
input.dimension.ARMA <- function(obj)
   {if (is.null(obj$C)) return(0)
    else   return(dim(obj$C)[3])
   }
input.dimension.TSestModel <- function(obj){input.dimension(obj$data)}

output.dimension <- function(x, ...)UseMethod("output.dimension")
output.dimension.SS <- function(obj){dim(obj$H)[1] }
output.dimension.ARMA <- function(obj){dim(obj$A)[2] }
output.dimension.TSestModel <- function(obj){output.dimension(obj$data)}




check.consistent.dimensions <- function(model,data)
   { # check data & model dimensions
    UseMethod("check.consistent.dimensions")
   }
check.consistent.dimensions.TSdata <- function(d,m)
   {check.consistent.dimensions(m,d)}
check.consistent.dimensions.TSestModel <- function(model, data)
   {if(missing(data)) data <- model$data
    check.consistent.dimensions(model$model,data)
   }

check.consistent.dimensions.SS <- function(model,data=NULL)
 {m <-dim(model$G)[2]
  n <-dim(model$F)[1]
  p <-dim(model$H)[1]
  if (n!= dim(model$F)[2]) stop("Model F matrix must be square.")
  if (n!= dim(model$H)[2])
      stop("Model H matrix have second dimension consistent with matrix F.")
  if (!is.null(model$G)) if(n!= dim(model$G)[1])
      stop("Model G matrix have first dimension consistent with matrix F.")
  if (!is.null(model$K)) if(n!= dim(model$K)[1])
      stop("Model K matrix have first dimension consistent with matrix F.")
  if (!is.null(model$K)) if(p!= dim(model$K)[2])
      stop("Model K matrix have second dimension consistent with matrix H.")

  if (!is.null(data))
   {if(dim(model$H)[1] != output.dimension(data))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(model$G))
      {if(0 != input.dimension(data))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(model$G)[2] != input.dimension(data))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(T)
 }

check.consistent.dimensions.ARMA <- function(model,data=NULL)
 {p <-dim(model$A)[2]
  if (p!= dim(model$A)[3]) stop("Model A array dim 2 and 3 should be equal.")
  if (p!= dim(model$B)[2]) stop("Model B array dim inconsistent with array A.")
  if (p!= dim(model$B)[3]) stop("Model B array dim inconsistent with array A.")
  if (!is.null(model$C))
      if (p!= dim(model$C)[2]) stop("Model C array dim inconsistent with array A.")

  if (!is.null(data))
   {if(dim(model$A)[2] != output.dimension(data))
       stop("Model and data output dimensions do not correspond.\n")
    if(is.null(model$C))
      {if(0 != input.dimension(data))
        stop("Model and data input dimensions do not correspond.\n")
      }
    else
      {if(dim(model$C)[3] != input.dimension(data))
         stop("Model and data input dimensions do not correspond.\n")
      }
   }
  return(T)
 }

check.consistent.dimensions.default <- function(model, ...)
   {stop(paste("No method for object of class ", dseclass(model), "\n"))}





TSestModel <- function(model) { model } # return everything





TSmodel <- function(model, ...) UseMethod("TSmodel")
TSmodel.TSmodel <- function(model){model}
#TSmodel.TSestModel <- function(emodel){emodel$model}

TSmodel.TSestModel <- function(model)
  {# Return a TSmodel object but also retains convergence info (if not null).
   model$model$converged <- model$converged
   model$model
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
   dseclass(model) <- c("ARMA","TSmodel")
   if(!is.null(names)) series.names(model) <- names
   else
     {if(!is.null( input.names))  input.series.names(model) <- input.names
      if(!is.null(output.names)) output.series.names(model) <- output.names
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
   if      (!is.null(model$K)) dseclass(model) <- c("innov","SS","TSmodel" )
   else if (!is.null(model$Q)) dseclass(model) <- c( "non.innov","SS","TSmodel")
   else stop("specified stucture is not correct for SS model.")
   if(!is.null(names)) series.names(model) <- names
   else
     {if(!is.null( input.names))  input.series.names(model) <-  input.names
      if(!is.null(output.names)) output.series.names(model) <- output.names
     }
   check.consistent.dimensions(model, data=NULL)
   set.parameters(model)
  }


parms <- function(model)UseMethod("parms")
parms.TSmodel <- function(model)  { model$parms }
parms.TSestModel <- function(model)  { model$model$parms }

is.TSmodel <- function(obj){inherits(obj,"TSmodel")}
is.TSestModel <- function(obj){inherits(obj,"TSestModel")}
is.SS <- function(obj){inherits(obj,"SS")}
is.innov.SS <- function(obj){inherits(obj,"SS")& inherits(obj,"innov")}
is.non.innov.SS <- function(obj){inherits(obj,"SS")&inherits(obj,"non.innov")}
is.ARMA <- function(obj){inherits(obj,"ARMA")}
is.TSm.or.em <- function(obj)
   {inherits(obj,"TSestModel") | inherits(obj,"TSmodel")}

TS.error.exit <- function(clss="TSmodel or TSestModel")
{ m <- sys.calls()
  mi <- m[[2]]
  m <- m[[length(m)-1]]
 stop(paste("Argument of class ",clss," required in call " 
   ,m[1],"(",m[2],") from initial call ", mi[1],"(",mi[2],")"))
}

# example:  
#  zot <- function(model) {  if(!is.TSm.or.em(model)) TS.error.exit()
#   2+2
#  }


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
   plist$parms <- c(plist$parms,A[indicate])          # parameters
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
    plist <- locateSS(model$F,model$const.F,"f",n,n,
                 list(parms=NULL,location=NULL,i=NULL,j=NULL,
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

    list.add(model, names(plist) ) <- plist
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
   plist$parms <- c(plist$parms,A[indicate])          # parameters
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
                     list(parms=NULL,location=NULL,i=NULL,j=NULL,
                         const=NULL,const.location=NULL,
                         const.i=NULL,const.j=NULL,l=NULL,const.l=NULL))
       plist <- locateARMA(model$B,"B",p,p,b,plist)
       if(!is.null(cc)) plist <- locateARMA(model$C,"C",p,m,cc,plist)
       if(!is.null(model$TREND)) 
            plist <- locateARMA(model$TREND,"t",p,m,cc,plist)

       list.add(model, names(plist) ) <- plist
       model
} #end set.parameters.ARMA


set.arrays <- function(model, parms=NULL)  
 { # complete representaion info. based on parameter info.  
  UseMethod("set.arrays")
 }
   
set.arrays.TSestModel <- function(model, parms=NULL)  
 {set.arrays(model$model, parms=parms) }
    
set.arrays.SS    <- function  (model, parms=NULL) 
{ # N.B. Dimension and class (innov/ non.innov) info. is assumed accurate
    if (is.null(parms)) parms   <-model$parms
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
    if (length(parms)>0) 
       {i <- a.pos == "f"
        f[cbind(i.pos[i],j.pos[i])] <- parms[i]
        if(!is.null(m)) 
          {i <- a.pos == "G"
           G[cbind(i.pos[i],j.pos[i])] <- parms[i]
          }
        i <- a.pos == "H"
        H[cbind(i.pos[i],j.pos[i])] <- parms[i]
        i <- a.pos == "K"
        K[cbind(i.pos[i],j.pos[i])] <- parms[i]
        i <- a.pos == "Q"
        Q[cbind(i.pos[i],j.pos[i])] <- parms[i]
        i <- a.pos == "R"
        R[cbind(i.pos[i],j.pos[i])] <- parms[i]
        i <- a.pos == "z"
        z[i.pos[i]] <- parms[i]
        i <- a.pos == "P"
        P[cbind(i.pos[i],j.pos[i])] <- parms[i]
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

set.arrays.ARMA     <- function  (model, parms=NULL) { 
# N.B. Dimension and class info. is assumed accurate
       if (is.null(parms)) parms   <-model$parms
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
       if (length(parms)>0) 
          {i <- a.pos == "A"
           A[cbind(l.pos[i],i.pos[i],j.pos[i])] <- parms[i]
           i <- a.pos == "B"
           B[cbind(l.pos[i],i.pos[i],j.pos[i])] <- parms[i]
           if(!is.null(m))
             {i <- a.pos == "C"
              C[cbind(l.pos[i],i.pos[i],j.pos[i])] <- parms[i]
             }
           i <- a.pos == "t"
           TREND[i.pos[i]] <- parms[i]
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

#######################################################################

#                    end

#######################################################################


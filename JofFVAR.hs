#!/bin/csh
#   2000/03/21 14:42:14 
#   NB   check load.from and  source("1dec93.dat") gives egJofF.1dec93.data 


# This is set up as a script file rather than an S function because the 
#  memory requirements are substantially less. (Splus 3.1 and 3.2 do not
#  release memory until the top level expression completes.)
#  NB dollar signs are changed to \$  for correct interpretation by the shell. If commands
#     are run interactively they should be changed back.
# The plots are not done by default. It is usually more practical to generate 
#   the objects with this script and then do the plots interactively.
#   See the end of this file for instructions to produce plots.

Splus <<eofS
 cat("This is sequence of instructions for producing results published in\n")
 cat("    Combining VAR Estimation and State Space Model Reduction\n")
 cat("       for Simple Good Predictions")
 cat("           Paul D. Gilbert\n")
 cat("  Journal of Forecasting, Special Issue on VAR Modelling and Forecasting\n")
 cat("           (forthcoming as of 1994)\n")
 cat(" Estimate of total time required is ( / / ) if mle is not run.\n")
 cat(" Timing estimates are based on (Sparc1+ / Sparc10/ Sparc1000)\n")
 cat(" Timings will depend substantially on the system load.\n")
 cat("\n")
 cat("NB. The results produced by these instuctions may vary slightly\n")
 cat("    from the published results. The random number generator had\n")
 cat("    a small change from Splus 3.1 to 3.2.\n")
 cat("\n")

 do.mle   <- F 
 do.plots <- F 

 #  data for section 6
 # source(".../1dec93.dat") # egJofF.1dec93.data  is loaded with DSE

 # check the random number generator (it changed from Splus 3.1 to Splus 3.2)
   assign(".Random.seed",c(37,60,43,61,39,3,4,38,24,48,49,1), where=1)
   z <- sum(abs(rnorm(100,sd=1)))
   orig <- 82.643575159296489119  #in Splus3.1 under Sun4
   if (z != orig) 
     {if (abs(z - orig) < 1e-15)
         {cat("Warning: The random number generator has small differences from the\n") 
          cat("version used to produce the published results.\n")
         }
      else if (abs(z - orig) < 1e-7)
         {cat("Warning: The random number generator has differences from the\n") 
          cat("version used to produce the published results.\n")
          }
      else 
         {cat("Warning: The random number generator has large differences from the\n") 
          cat("version used to produce the published results.\n")
          }
      cat("\n")
     }
 rm(z, orig)



 mod1 <- ARMA(A=array(c(1,-.25,-.05), c(3,1,1)), B=array(1,c(1,1,1)))
 mod2 <- ARMA(A=array(c(1,-.8, -.2 ), c(3,1,1)), B=array(1,c(1,1,1)))
 mod3 <- ARMA(A=array(c( 
       1.00,-0.06,0.15,-0.03,0.00,0.02,0.03,-0.02,0.00,-0.02,-0.03,-0.02,
       0.00,-0.07,-0.05,0.12,1.00,0.20,-0.03,-0.11,0.00,-0.07,-0.03,0.08,
       0.00,-0.40,-0.05,-0.66,0.00,0.00,0.17,-0.18,1.00,-0.11,-0.24,-0.09 ) ,c(4,3,3)), 
      B=array(diag(1,3),c(1,3,3)))

  cat("J of F VAR paper section 2 generating object 1 (10 min/ 7 min /5 min)\n")
  cat("     started ",date.parsed(), "...\n")

  e.ls.mod1 <- eval.estimation( mod1, replications=100, 
         seed=c(13,44,1,25,56,0,6,33,22,13,13,0),
         simulation.args=list(sampleT=100, sd=1), 
         estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
          criterion="TSmodel")

  cat("J of F VAR paper section 2 generating object 2 (10 min/ 7 min / 5 min)\n")
  cat("     started ",date.parsed(), "...\n")
  e.ls.mod2 <- eval.estimation( mod2, replications=100, 
         seed=c(13,43,7,57,62,3,30,29,24,54,47,2),
         simulation.args=list(sampleT=100, sd=1), 
         estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
          criterion="TSmodel")

  cat("J of F VAR paper section 2 generating object 3 (15 min/ 8 min / 10 min)\n")
  cat("     started ",date.parsed(), "...\n")
  e.mod3.VAR.lag12<-eval.estimation( mod3, replications=100, 
         seed=c(37,17,6,24,47,2,62,7,62,2,21,3),
         simulation.args=list(sd=1), 
         estimation="est.VARX.ls", estimation.args=list(max.lag=12), 
          criterion="TSmodel")

  cat("J of F VAR paper section 2 generating object 4 (20 min/ 8 min /10 min)\n")
  cat("     started ",date.parsed(), "...\n")
  e.mod3.VAR.lag6<-eval.estimation( mod3, replications=100, 
         seed=c(53,41,26,39,10,1,19,25,56,32,28,3),
         simulation.args=list(sd=1), 
         estimation="est.VARX.ls", estimation.args=list(max.lag=6), 
         criterion="TSmodel")

  cat("J of F VAR paper section 2 generating object 5 (10 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.mod3.VAR.lag12 <- forecast.cov.estimators.wrt.true(mod3,
     seed=c(37,17,6,24,47,2,62,7,62,2,21,3), 
     estimation.methods=list(est.VARX.ls=list(max.lag=12)),
     est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 2 generating object 6 (5 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), " ...\n")
  pc.mod3.VAR.lag6 <- forecast.cov.estimators.wrt.true(mod3,
     seed=c( 53,41,26,39,10,1,19,25,56,32,28,3), 
     estimation.methods=list(est.VARX.ls=list(max.lag=6)),
     est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 3 generating object 1 ( / 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.rd.ls.3lag <- forecast.cov.reductions.wrt.true(mod3,
     seed=c(29,55,47,18,33,1,15,15,34,46,13,2),
     estimation.methods=list(est.VARX.ls=list(max.lag=3)),
     est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 3 generating object 2 (4.5 hrs/ 25 min /25 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.rd.ls.12lag <- forecast.cov.reductions.wrt.true(mod3,
     seed= c(53,33,11,11,54,3,54,15,33,28,9,2),
     estimation.methods=list(est.VARX.ls=list(max.lag=12)),
      est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 3 generating object 3 (10 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.ewt7ls.3.12lag <- forecast.cov.estimators.wrt.true(mod3, 
     seed=c(13,61,61,38,23,1,63,44,34,19,59,2),
     estimation.methods=list(
           est.VARX.ls=list(max.lag=3, lag.weight=.7),
           est.VARX.ls=list(max.lag=12, lag.weight=.7)),
     est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 3 generating object 4 (2 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.ewt7ls.mod2.12lag <- forecast.cov.estimators.wrt.true(mod2,
      seed=c(37,60,43,61,39,3,4,38,24,48,49,1),
      estimation.methods=list(est.VARX.ls=list(max.lag=12, lag.weight=.7)),
      est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 3 generating object 5 (3 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.ewt9ls.mod2.12lag <- forecast.cov.estimators.wrt.true(mod2,
      seed=c(37,60,43,61,39,3,4,38,24,48,49,1),
      estimation.methods=list(est.VARX.ls=list(max.lag=12, lag.weight=.9)),
      est.replications=2, pred.replications=10, Spawn=T) 

  cat("J of F VAR paper section 3 generating object 6 (15 min/ 5 min /5 min)\n")
  cat("     started ",date.parsed(), "...\n")
  e.ewt7ls.mod2.12lag <- eval.estimation( mod2, replications=100,
     seed=c(29,55,47,18,33,1,15,15,34,46,13,2),
     estimation="est.VARX.ls",
     estimation.args=list(max.lag=12, lag.weight=.7),
     criterion="TSmodel")

  cat("J of F VAR paper section 4 generating object 1 ( / 20 min /20 min)\n")
  cat("     started ",date.parsed(), "...\n")

    pc.rd.ewt7ls.12lag <- forecast.cov.reductions.wrt.true(mod3,
     seed=c(29,16,40,58,14,2,41,2,38,24,56,0),
     estimation.methods=list(est.VARX.ls=list(max.lag=12, lag.weight=.7)),
     est.replications=2, pred.replications=10, Spawn=T, criteria="taic")

  cat("J of F VAR paper section 4 generating object 2 (6 min/ 2 min /2 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.ls.mod3.known <-  forecast.cov.estimators.wrt.true(mod3,
     seed=c(29,16,40,58,14,2,41,2,38,24,56,0),
     estimation.methods=list(
           est.VARX.ls=list(max.lag=3),
           est.VARX.ls=list(max.lag=3, lag.weight=.7)),
     est.replications=2, pred.replications=10, Spawn=T)

   mod4 <- mod3
   mod4\$B <- array(0,c(2,3,3))
   mod4\$B[1,,] <- diag(1,3)
   mod4\$B[2,,] <- diag(0.9,3)
   mod4 <- set.parameters(mod4)

   mod5 <- mod4
   zl<- array(0,c(2,3,3))
   zl[1,,] <- diag(1,3)
   zl[2,1,1] <- -1
   mod5\$A <- polyprod(zl, mod4\$A)
   mod5 <- set.parameters(mod5)
   all.ok <- T

  cat("J of F VAR paper section 5 generating object 1 ( / /days)\n")
  cat("     started ",date.parsed(), "...\n")
  if (do.mle)
    {  cat("(\nThis object has only been computed on a Sparc10\n")
       pc.mle.mod4.known <-  forecast.cov.estimators.wrt.true(mod4,  
       seed=c(17,13,42,31,5,3,10,52,17,44,54,3),
       estimation.methods=list(
           est.max.like=list(mod4,max.iter=50, algorithm="dfpMin",
                             line.search="brent", dfunc=numerical.grad)), 
       est.replications=2, pred.replications=10, Spawn=T)

     # the following "converges", but to an inferior value than above.
     # first estimate iteration  7200+ value  411.893025669975 (conv.)
     # second         iteration  9100+ value  410.890331230053 (conv.)
       nls.pc.mle.mod4.known <-  forecast.cov.estimators.wrt.true(mod4,  
       seed=c(17,13,42,31,5,3,10,52,17,44,54,3),
       estimation.methods=list(
           est.max.like=list(mod4,max.iter=10000, algorithm="nlsimplex")), 
       est.replications=2, pred.replications=10, Spawn=T)
    }
#  else cat("skipped\n")

  cat("J of F VAR paper section 5 generating object 2 ( / 30 min /30 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.bb4.mod4       <-  forecast.cov.estimators.wrt.true(mod4,  
    seed=c(17,13,42,31,5,3,10,52,17,44,54,3), 
    estimation.methods=list(
          est.black.box4=list(verbose=F, criterion="taic", max.lag=12)), 
    est.replications=2, pred.replications=10, Spawn=T)

  cat("J of F VAR paper section 5 generating object 3  ( / 30 min /30 min)\n")
  cat("     started ",date.parsed(), "...\n")
  pc.bb4.mod5       <-  forecast.cov.estimators.wrt.true(mod5,  
    seed=c(17,13,42,31,5,3,10,52,17,44,54,3), 
    estimation.methods=list(
           est.black.box4=list(verbose=F, criterion="taic", max.lag=12)), 
    est.replications=2,
    pred.replications=10, Spawn=T)

# This assumes there is a version of required compiled code, compiled with a 
# larger state dimension, in the search path were DSE is found.
load.from <-"/home/res/gilp/Sroutines"
#load.from <- search()[grep("b*/DSE/.Data",search())]
#load.from <- substring(load.from,1,nchar(load.from)-6)
load.DSE.large.fortran(from=load.from)
 if(exists(".First", where=1)) orig.First <- .First
.First <- function()
{# this is necessary to do the following to 
 # with max.lag=12. max.lag=10 should be possible without this.
 if(exists("orig.First", where=1)) orig.First()
 invisible(load.DSE.large.fortran(from=load.from))
}

   cat("J of F VAR paper section 6 generating object 1 (/ 45 min /45 min)\n")
   cat("     started ",date.parsed(), "...\n")
   bb4.vs.ls.cov <-out.of.sample.forecast.cov.estimators.wrt.data(egJofF.1dec93.data, 
       zero=T, trend=T, estimation.sample=.75, 
       estimation.methods = list(
           est.black.box4=list(estimation="est.VARX.ls", max.lag=12, verbose=F),
                    # use verbose=T above to get info. criteria printout 
           est.VARX.ls=list(max.lag=12),
           est.VARX.ls=list(max.lag=6)))  

 if(exists("orig.First", where=1))
   {.First <- orig.First
    rm(orig.First)
   }

   cat("J of F VAR paper section 6 generating object 2 ( / 20 min /23 min)\n")
   cat("     started ",date.parsed(), "...\n")
   bb4.vs.ls.proj <-horizon.forecasts(bb4.vs.ls.cov,horizons=c(1,3,6),
                        discard.before=1)
   egJofF.1dec93.data.subset<- egJofF.1dec93.data
   egJofF.1dec93.data.subset\$output <- 
           egJofF.1dec93.data.subset\$output[,c(1,2,6,7)]
   tframe(egJofF.1dec93.data.subset\$output) <- 
           tframe(egJofF.1dec93.data\$output)

   subset.bb4.vs.ls.cov<-out.of.sample.forecast.cov.estimators.wrt.data( 
      egJofF.1dec93.data.subset, zero=T, trend=T, estimation.sample=.75, 
      estimation.methods = list( 
         est.black.box4=list(  estimation="est.VARX.ls", max.lag=12, verbose=F),
         est.VARX.ls=list(max.lag=12),
         est.VARX.ls=list(max.lag=6) ))

  cat("J of F VAR paper - object generation finished ",date.parsed(), "\n")

  if (do.plots) 
    {# it is necessary to have a graphic device for plotting, eg:
          #  openlook()
          #  sunview()
          #  motif()
          #  or use postscript to create graphics file:
          #     postscript(file="figx.ps",width=6,height=6,pointsize=10,
          #           onefile=F, print.it=F, append=F)
          #     tfplot(...)
          #     dev.off()

     par(mfcol=c(2,1))
     tfplot(parms(e.ls.mod1), bounds=F)    # line chart of cum. ave. est. parameters
     tfplot(parms(e.ls.mod2), bounds=F)    # line chart of cum. ave. est. parameters
     distribution(parms(e.ls.mod1), bandwidth=.2)  # plot of dist. of estimates
     distribution(parms(e.ls.mod2), bandwidth=.2)  # plot of dist. of estimates
     plot(roots(e.ls.mod1), complex.plane=F)   # line chart of cum. ave. est. roots
     plot(roots(e.ls.mod2), complex.plane=F)   # line chart of cum. ave. est. roots
     distribution(roots(e.ls.mod1), bandwidth=.2)  #Exhibit 1
     distribution(roots(e.ls.mod2), bandwidth=.1)  #Exhibit 2
     tfplot(total.forecast.cov(pc.mod3.VAR.lag12), select.trend=F)  #Exhibit 3
     tfplot(total.forecast.cov(pc.mod3.VAR.lag6 ), select.trend=F)  #Exhibit 4
     tfplot(pc.mod3.VAR.lag6, select.trend=F, cex=1)                #Exhibit 5

     tfplot(total.forecast.cov(pc.rd.ls.3lag), select.trend=F, 
                                      lty=c(1,2,rep(3,18)) )     #Exhibit 6
      #legend(6,4.1,c("true", "zero", "sum of squared error"), lty=1:3)
     tfplot(total.forecast.cov(pc.rd.ls.12lag), select.trend=F,
                                      lty=c(1,2,rep(3,36)) )     #Exhibit 7
     tfplot(total.forecast.cov(pc.ewt7ls.3.12lag),select.trend=F,
                                      lty=c(1,2,rep(3,18)) )     #Exhibit 8
     tfplot(pc.ewt7ls.mod2.12lag,select.trend=F)                   #Exhibit 9
     tfplot(pc.ewt9ls.mod2.12lag,select.trend=F)                   #Exhibit 10
     distribution(roots(e.ewt7ls.mod2.12lag,
          criterion.args = list(randomize = T)), select= 1)      #Exhibit 11

     tfplot(total.forecast.cov(pc.rd.ewt7ls.12lag), select.trend=F,
                                      lty=c(1,2,1,rep(3,36)) )   #Exhibit 12
     print(pc.rd.ewt7ls.12lag\$information.criteria[[2]], digits=3)  #Exhibit 13

     tfplot(total.forecast.cov(
         combine(pc.ls.mod3.known, pc.rd.ewt7ls.12lag)),
         select.trend=F, lty=c(1,2,1,1,rep(3,36)) )              #Exhibit 14

     tfplot(total.forecast.cov(combine(pc.bb4.mod4,pc.mle.mod4.known)), 
            select.trend=F, lty=c(1,2,3,4) )   #Exhibit 15

     tfplot(total.forecast.cov(pc.bb4.mod5), select.zero=F, select.trend=F,
                                      lty=c(1,2,rep(3,18)) )   #Exhibit 16

#     selected series forecasts starting in June 1989 for  
#          horizons of 1, 3 and 6  months -  bft model 
     tfplot(bb4.vs.ls.proj[[1]], select.series=c(1,2,6,7), start.=c(1989,6)) 

#     selected series forecasts starting in June 1989 for  
#          horizons of 1, 3 and 6  months -  12 lag VAR 
     tfplot(bb4.vs.ls.proj[[2]], select.series=c(1,2,6,7), start.=c(1989,6)) 

#     selected series forecasts starting in June 1989 for  
#          horizons of 1, 3 and 6  months -  6 lag VAR 
     tfplot(bb4.vs.ls.proj[[3]], select.series=c(1,2,6,7), start.=c(1989,6)) 

     tfplot(bb4.vs.ls.cov, select.series=c(1,2,6,7),select.zero=F, 
          select.trend=F, lty=c(1,2,3), cex=1, mar=c(5,8,4,2)) #Exhibit 17

     tfplot(bb4.vs.ls.cov, select.series=c(1,2,6,7),select.cov=c(1,3),
        select.zero=F, select.trend=T, cex=1, mar=c(5,8,4,2))  #Exhibit 18

# plots of residual diagnostics
# bft model:
# check.residuals(l(bb4.vs.ls.cov\$multi.model[[1]],bb4.vs.ls.cov\$data))
# VAR model:
# check.residuals(l(bb4.vs.ls.cov\$multi.model[[2]],bb4.vs.ls.cov\$data))

     tfplot(subset.bb4.vs.ls.cov,lty=c(1,2,3),select.cov=c(1,3),
       select.zero=F,select.trend=T, cex=1, mar=c(5,8,4,2))    #Exhibit 19

# plots of residual diagnostics
# small bft model:
# check.residuals(l(subset.bb4.vs.ls.cov\$multi.model[[1]],subset.bb4.vs.ls.cov\$data))

# small VAR model:
# check.residuals(l(subset.bb4.vs.ls.cov\$multi.model[[2]],subset.bb4.vs.ls.cov\$data))


    }

q()
eofS
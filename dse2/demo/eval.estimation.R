#######  examples from User's Guide section 

  require("dse2")
 
#  make a "true" model to evaluate estimation algorithms

  mod1 <- ARMA(A=array(c(1,-.25,-.05), c(3,1,1)), B=array(1,c(1,1,1)))
  mod2 <- ARMA(A=array(c(1,-.8, -.2 ), c(3,1,1)), B=array(1,c(1,1,1)))
  mod3 <- ARMA(
 	A=array(c( 
 	1.00,-0.06,0.15,-0.03,0.00,0.02,0.03,-0.02,0.00,-0.02,-0.03,-0.02,
	0.00,-0.07,-0.05,0.12,1.00,0.20,-0.03,-0.11,0.00,-0.07,-0.03,0.08,
 	0.00,-0.40,-0.05,-0.66,0.00,0.00,0.17,-0.18,1.00,-0.11,-0.24,-0.09 )
		,c(4,3,3)), 
 	B=array(diag(1,3),c(1,3,3)))
  e.ls.mod1 <- EstEval( mod1, replications=100, 
 	simulation.args=list(sampleT=100, sd=1), 
 	estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
 	criterion="TSmodel", quiet=T)

     tfplot(coef(e.ls.mod1))
     tfplot(coef(e.ls.mod1), cum=F, bounds=F) 
     distribution(coef(e.ls.mod1), bandwidth=.2)


 #  this may take awhile
 
  e.ls.mod2 <- EstEval( mod2, replications=100, 
                     simulation.args=list(sampleT=100, sd=1), 
                     estimation="est.VARX.ls", estimation.args=list(max.lag=2), 
                     criterion="TSmodel", quiet=T)

     
     tfplot(coef(e.ls.mod2)) 
  Sys.sleep(3)
     tfplot(coef(e.ls.mod2), cum=F, bounds=F) 
  Sys.sleep(3)
     distribution(coef(e.ls.mod2), bandwidth=.2)
    
  
  e.ls.mod1.roots <- roots(e.ls.mod1)
  
     plot(e.ls.mod1.roots) 
  Sys.sleep(3)
     plot(e.ls.mod1.roots, complex.plane=F)
     
     plot(roots(e.ls.mod2), complex.plane=F) 
  Sys.sleep(3)
     distribution(e.ls.mod1.roots, bandwidth=.2) 
  Sys.sleep(3)
     distribution(roots(e.ls.mod2), bandwidth=.1) 
  

  pc <- forecastCov.estimators.wrt.true(mod3,
 	estimation.methods=list(est.VARX.ls=list(max.lag=6)),
 	est.replications=2, pred.replications=10, quiet=T)

  tfplot(pc)
 
  pc.rd <- forecastCov.reductions.wrt.true(mod3,
 	estimation.methods=list(est.VARX.ls=list(max.lag=3)),
 	est.replications=2, pred.replications=10, quiet=T)

  tfplot(pc.rd)


  data("eg1.DSE.data", package="dse1")

  z <-out.of.sample.forecastCov.estimators.wrt.data(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.VARX.ar=list(warn=F), 
 		est.VARX.ls=list(warn=F)), 
 	trend=T, zero=T)
  tfplot(z)

  zz <-out.of.sample.forecastCov.estimators.wrt.data(trim.na(eg1.DSE.data),
 	estimation.sample=.5,
 	estimation.methods = list(
 		est.black.box4=list(max.lag=3, verbose=F, warn=F),
		est.VARX.ls=list(max.lag=3, warn=F)), 
	trend=T, zero=T)

  tfplot(zz)
  
  
  zf<-horizonForecasts(TSmodel(zz, select=1),TSdata(zz), horizons=c(1,3,6))

  tfplot(zf)

#######  end of examples from User's Guide section 

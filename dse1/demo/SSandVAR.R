
require("dse1")

 data("eg1.DSE.data", package = "dse1")
 data("eg1.DSE.data.diff", package = "dse1") 

 cat("Retrieve data from file eg1.dat and build a TSdata object.\n")

#eg1.dat <- t(matrix(scan("eg1.dat", sep="")),5, 364))[, 2:5] 
#eg1.dat <- TSdata(input  = eg1.dat[,1,drop = F], 
#                  output = eg1.dat[, 2:4, drop = F])
#eg1.dat <-tframed(eg1.dat, list(start=c(1961,3), frequency=12))

eg1.dat <-tframed(eg1.DSE.data, list(start=c(1961,3), frequency=12))

seriesNamesInput(eg1.dat) <- "R90"
seriesNamesOutput(eg1.dat) <- c("M1","GDPl2", "CPI")

  cat("This should be the same as eg1.DSE.data\n")
if (test.equal(eg1.dat, eg1.DSE.data))
  cat("Data compares ok.\n") else
  cat("Data does not compare.\n")


  cat("Note: The first part of this demo uses the total sample. BoC Working\n")
  cat("    Paper 93-4 was estimated with a sub-sample and used est.VARX.ar!\n")


#  VARmodel <- est.VARX.ar(eg1.DSE.data.diff, re.add.means=F)
  VARmodel <- est.VARX.ls(eg1.DSE.data.diff, re.add.means=F)
  summary(VARmodel)
  roots(VARmodel)
  stability(VARmodel)
  tfplot(VARmodel)

#  cat("Use a TSmodel and evaluate it using a \n")
#  cat("different data set (the first 240 observations from eg1.DSE.data.diff).\n")
# summary(l(VARmodel, eg1.DSE.data.diff, sampleT=240, predictT=240))


  cat("Use a TSmodel and evaluate it using a different data set.\n")
  summary(l(VARmodel, eg1.DSE.data.diff, sampleT=363, predictT=363))

  cat("       to.SS and evaluate with the same data set...\n")
  SSmodel <- l(to.SS(VARmodel), TSdata(VARmodel)) 
  abs(VARmodel$estimates$like[1] - SSmodel$estimates$like[1])
 

  cat("     to.ARMA and evaluate with the same data set...\n")
  ARMAmodel <- l(to.ARMA(SSmodel), TSdata(VARmodel))
  abs(VARmodel$estimates$like[1] - ARMAmodel$estimates$like[1])
  

  cat(" McMillan.degree...\n")
  McMillan.degree(VARmodel, verbose=F)$distinct
   

  cat(" Working Paper 93-4 comparisons:\n")
  cat("      VAR model likelihood...\n")
  sub.sample <- tfwindow(eg1.DSE.data.diff,end=c(1981,2))
  VARmodel <- est.VARX.ar(sub.sample, re.add.means=F)
  summary(VARmodel) 

  cat("   VAR model \n")
  roots(VARmodel, by.poly=T)

  cat("   State Space  model\n")
  SS1.model <- l(balance.Mittnik(to.SS(VARmodel), n=9), sub.sample)
  summary(SS1.model)
   
  roots(SS1.model)

  cat("   ARMA model \n")
  ARMA.model<- l(to.ARMA(SS1.model), sub.sample)
  summary(ARMA.model) 

  cat("     ARMA model roots  (calculated two ways) ...\n")
  roots(ARMA.model, fuzz=1e-4)
  roots(ARMA.model, fuzz=1e-4, by.poly=T)
  
  stability(ARMA.model, fuzz=1e-4)


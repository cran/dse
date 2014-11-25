#######  examples from User's Guide section 6

  require("dse1")
  
  if (is.R()) data("eg1.DSE.data", package = "dse1") else 
 if (is.S()) 
   {source(paste(DSE.HOME, "/data/eg1.DSE.data.R", sep=""))
    class(eg1.DSE.data$output) <- class(eg1.DSE.data$input) <- NULL
    }

  cat("Estimation...\n")

  model.eg1.ls <- est.VARX.ls(trim.na(eg1.DSE.data), warn=F)
  subsample.data <- tfwindow(eg1.DSE.data,start=c(1972,1),end=c(1992,12))

  summary(model.eg1.ls)
  model.eg1.ls # or print(model.eg1.ls)
  
  tfplot(model.eg1.ls)
  tfplot(model.eg1.ls, start.=c(1990,1))

  check.residuals(model.eg1.ls, plot.=F, pac=T)


  model.eg1.ss <- est.SS.from.VARX(trim.na(eg1.DSE.data)) 
  summary(model.eg1.ss)
  model.eg1.ss # or print(model.eg1.ss)
 
  information.tests(model.eg1.ls, model.eg1.ss)
 

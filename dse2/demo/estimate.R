#######  examples from User's Guide section 6

  require("dse2")
  
  data("eg1.DSE.data", package = "dse1")

##### Estimation

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
 
  data("egJofF.1dec93.data", package="dse1")
  #  select a subset of the data
  eg4.DSE.data<- egJofF.1dec93.data
  output.data(eg4.DSE.data) <- output.data(eg4.DSE.data, series=c(1,2,6,7))

  model.eg4.bb <- est.black.box(trim.na(eg4.DSE.data), max.lag=3, verbose=F) 

  tfplot(model.eg4.bb)



#######  end of examples from User's Guide section 6

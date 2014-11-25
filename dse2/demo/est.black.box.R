
  require("dse2")
  
 if (is.R()) data("egJofF.1dec93.data", package="dse1")else 
 if (is.S()) 
   {source(paste(DSE.HOME, "/data/egJofF.1dec93.data.R", sep=""))
    class(egJofF.1dec93.data$output) <- class(egJofF.1dec93.data$input) <- NULL
    }

  cat("select a subset of the data\n")
  eg4.DSE.data<- egJofF.1dec93.data
  output.data(eg4.DSE.data) <- output.data(eg4.DSE.data, series=c(1,2,6,7))

  cat("Estimation...\n")
  model.eg4.bb <- est.black.box(trim.na(eg4.DSE.data), max.lag=3, verbose=F) 

  tfplot(model.eg4.bb)
  summary(model.eg4.bb)

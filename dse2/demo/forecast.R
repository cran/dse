#######  examples from User's Guide section 7

  require("dse2")
  
  data("egJofF.1dec93.data", package="dse1")
  eg4.DSE.data<- egJofF.1dec93.data
  output.data(eg4.DSE.data) <- output.data(eg4.DSE.data, series=c(1,2,6,7))

#  Forecasting

  eg4.DSE.model <- est.VARX.ls(eg4.DSE.data)
 
  new.data <- TSdata(
              input= ts(rbind(input.data(eg4.DSE.data), matrix(.1,10,1)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)),    
              output=ts(rbind(output.data(eg4.DSE.data),matrix(.3,5,4)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)))
  series.names(new.data) <- series.names(eg4.DSE.data)

  z  <- l(TSmodel(eg4.DSE.model), trim.na(new.data)) 

  zz <- forecast(TSmodel(eg4.DSE.model), new.data)
  z <-  forecast(TSmodel(eg4.DSE.model), trim.na(new.data), 
		conditioning.inputs=input.data(new.data))
  tfplot(zz, start.=c(1990,6))


  z <- forecast(eg4.DSE.model, conditioning.inputs.forecasts=matrix(.5,6,1)) 
  summary(z)
  tfplot(z)
  Sys.sleep(3)
  tfplot(z, start.=c(1990,1))
  

  z <- l(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  
  z <- featherForecasts(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  
  zz <-featherForecasts(TSmodel(eg4.DSE.model), new.data,
                          from.periods =c(20,50,60,70,80), horizon=150)
  tfplot(zz)


  z <- horizonForecasts(TSmodel(eg4.DSE.model), new.data, horizons=c(1,3,6))
  tfplot(z)


  fc1 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data)
  
  tfplot(fc1)
  Sys.sleep(3)
  tfplot(forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, horizons= 1:4)) 
 
  fc2 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, zero=T, trend=T)
  tfplot(fc2)  # may want mar=c(2.1,4.1,4.1,2.1)

#######  end of examples from User's Guide section 7

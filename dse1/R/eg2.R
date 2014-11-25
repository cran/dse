
example.VAR.SVD <- function(d,d.all, d.all.raw)
 {V.1 <- est.VARX.ar(d) # estimates a VAR model using the truncated sample.

  l.V.1 <-l(V.1, d) # calculates the likelihood, one step ahead predictions, etc., and puts them in the variable l.V.1.

  cat("Likelihood and components for VAR model\n")
  # prints the likelihood value (with a breakdown for the 3 terms of the likelihood function).
  print(l.V.1$estimates$like, digits=16)

  #cat("Likelihood and components for VAR model\n")
  #l.V.1$estimates$like # also prints the value but not as many digits.
  # calculate the likelihood, one step ahead predictions, etc., based 
  #   on the full sample, and puts them in the variable o.V.1.
  o.V.1 <-l(V.1, d.all) 

  # convert the VAR model to a state space model balanced by Mittnik's technique.
  SS.V.1 <- to.SS(V.1) 

  # calculate the likelihood, one step ahead predictions, etc., based on 
  #   the truncated sample, and puts them in the variable l.SS.V.1.
  l.SS.V.1 <-l(SS.V.1, d) 

  cat("Likelihood and components for state space model\n")
  # print the likelihood value (with a breakdown for the 3 terms of 
  #     the likelihood function).
  print(l.SS.V.1$estimates$like,digits=16) 

  cat("Maximum difference in one-step predictions of VAR and state space model ")
  # calculate the difference of the absolute values of the predictions of 
  #      the two models.
  cat(max(abs(l.V.1$estimates$pred - l.SS.V.1$estimates$pred))) 
  cat("\n")

  cat("Exhibit 2. Mittnik reduction from VAR model: \n")
  M5.SS.V.1 <- reduction.Mittnik(SS.V.1,data=d, criterion="taic")  
   #If criterion is not specified the program prompts for a state dimension and 
   #returns that model. Results is put in the variable M5.SS.V.1.

  cat("Exhibit 3. Mittnik estimation lag=3: \n")
  M12.shift3 <- est.SS.Mittnik(d,max.lag=3, n=12)
  M12.shift3 <- reduction.Mittnik(M12.shift3, data=d, criterion="taic")  

  cat("Exhibit 4. Mittnik estimation lag=4: \n")
  M12.shift4 <- est.SS.Mittnik(d,max.lag=4, n=15)
  M12.shift4 <- reduction.Mittnik(M12.shift4, data=d, criterion="taic")  

  if( dev.cur() != 1 ) # eg.-OpenLook() or Suntools() 
     example.show.ytoy.cpi (o.V.1,d.all.raw,start=240)
  else  cat("Exhibit 8. graphic requires graphic device. \n")
  invisible()
}



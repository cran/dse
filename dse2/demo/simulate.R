#######  examples from User's Guide section 5

  require("dse2")

####  ARMA model

  ar   <- array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3), c(3,2,2))
  ma   <- array(c(1,.2,0,.1,0,0,1,.3), c(2,2,2))
  arma <- ARMA(A=ar, B=ma, C=NULL) # make an ARMA TSmodel object

  arma  # or print(arma)

  data.arma.sim <- simulate(arma)  # generate simulated data
  arma <- l(arma,data.arma.sim)    # evaluate the model with the simulated data
                                   # to get a TSestModel object
  summary(arma)
  roots(arma)
  stability(arma)

  tfplot(data.arma.sim)
  tfplot(arma)

####  State Space model

  f  <- array(c(.5,.3,.2,.4),c(2,2))
  h  <- array(c(1,0,0,1),c(2,2))
  k  <- array(c(.5,.3,.2,.4),c(2,2))
  ss <- SS(F=f,G=NULL,H=h,K=k)      # make an SS TSmodel object
  ss   # or print(ss)
  
  data.ss.sim <- simulate(ss)   # generate simulated data
  ss <- l(ss,data.ss.sim)       # evaluate the model with the simulated data
                                # to get a TSestModel object
  summary(ss)
  roots(ss)
  stability(ss)
  
  tfplot(ss)

####  Convert between State Space and ARMA models

  ss.from.arma <- l(to.SS(arma), data.arma.sim)
  arma.from.ss <- l(to.ARMA(ss), data.ss.sim)
  
  summary(ss.from.arma)
  summary(arma.from.ss)
  
  stability(arma)
  stability(ss.from.arma)
  
#######  end of examples from User's Guide section 5 

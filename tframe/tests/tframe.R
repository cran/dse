  require("tframe")

 Sys.info()

  #x11()
#  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
#             # width=6, height=8, pointsize=10,



tframe.function.tests <- function( verbose=T, synopsis=T)
{# A short set of tests of the tframe class methods. 

  all.ok <-  T
  if (synopsis & !verbose) cat("All tframe tests ...")
  if (verbose) cat("tframe test 1 ... ")
  tspvector <- tframed(1:100, list(start=c(1981,3), frequency=4))
  data <- matrix(rnorm(300),100,3)
  tframe(data) <- tframe(tspvector)
  ok <- is.tframed(data)
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("tframe test 2 ... ")
  ok <- test.equal(tframe(data), tframe(data))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 3 ... ")
  ok <- all(c(1981,3) == start(tspvector))
  ok <- ok & all(c(1981,3) == start(data))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 4 ... ")
  ok <- all(end(data) == end(tspvector))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 5 ... ")
  ok <- periods(data) == periods(tspvector)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 6 ... ")
  ok <- frequency(data) == frequency(tspvector)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 7 ... ")
  z <- tframed(data, list(start=c(1961,2), frequency=12) )
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 8 ... ")
  z <- data[10:90,]
  tframe(z) <- tfTruncate(tframe(data), start=10, end=90)
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 9 ... ")
  z <- tfTruncate(data, start=10, end=90)
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 10... ")
  data <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  z   <- tfwindow(data, start=c(1963,2))
  zz  <- data
  zz  <- tfwindow(zz, start=c(1963,2))
  zzz <- tfwindow(data, start=c(1963,2))
  tframe(zzz) <- tframe(z)
  zzz  <- tframed(zzz, tframe(zzz))
  zzzz <- tframed(matrix(rnorm(300),100,3), tframe(data))

  all.ok <- all.ok &  is.tframed(data)
  all.ok <- all.ok &  is.tframed(z)
  all.ok <- all.ok &  is.tframed(zz)
  all.ok <- all.ok &  is.tframed(zzz)
  all.ok <- all.ok &  is.tframed(zzzz)
  all.ok <- all.ok &  all(z==zz)
  all.ok <- all.ok &  all(z==zzz)
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 11... ")
  ok <- all( time(data) == time( tframed(data, tframe(data))))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 12... ")
  z <- tsmatrix(1:10, 11:20)
  ok <-  all(start(z) ==1) & all( z== matrix(1:20, 10,2)) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 13... ")
  data <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  z <- tfwindow(data, start=c(1963,2), end=c(1969,1))
  ok <-      all(start(data)== earliestStart(data, z))
  ok <- ok & all(    end(z) == earliestEnd  (data, z))
  ok <- ok & all(start(z)   == latestStart  (data, z))
  ok <- ok & all( end(data) == latestEnd   (data, z))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 14... ")
  data <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  z <- tfwindow(data, start=c(1963,2), end=c(1969,1))
  ok <- test.equal(data, splice(z, data))
  ok <- ok & test.equal(tframe(data), tframe(splice(z, data)))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if(dev.cur() == 1) postscript(file="tframeGraphicsTest.ps")

# plot(data)
  tfplot(data)
  seriesNames(data) <- c("newname 1", "newname 2", "newname 3")
  tfplot(data)

  if (synopsis) 
    {if (verbose) cat("All tframe tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n") }
     
  if (all.ok) invisible(T)  else stop("FAILED")
}


   tframe.function.tests(verbose=T) 

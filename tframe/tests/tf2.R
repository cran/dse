 require("ts")
 require("tframe")

 Sys.info()

t1<-ts(c(1,2,3,4,5),start=c(1991,1))
t2<-ts(c(2,3,4,5,6,7,8),start=c(1992,1))
t3<-ts(c(NA,2,3,4,5),start=c(1991,1))

latest.start(t1,t2,t3) # 1992 1 corresponding to the starting date of 
                       # the object which starts latest (t2)
ok <- all(c(1992,1) == latest.start(t1,t2,t3))
ok

latest.start(t1,t3)     # both start in 1991 1 (NAs count as data)
ok <- ok & all(c(1991,1) == latest.start(t1,t3))
ok

latest.start(tbind(t1,t2,t3)) # tbind gives a single object starting in 1991 1
ok <- ok & all(c(1991,1) == latest.start(tbind(t1,t2,t3)))
ok

latest.start(t2, tbind(t1,t2,t3))
ok <- ok & all(c(1992,1) == latest.start(t2, tbind(t1,t2,t3)))
ok

latest.start.index(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 2 == latest.start.index(t1,t2,t3)
ok

latest.end.index(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 2 == latest.end.index(t1,t2,t3)
ok

earliest.end.index(t1,t2,t3)  # position of t2 in the argument list
ok <- ok & 1 == earliest.end.index(t1,t2,t3)
ok

earliest.start(t1,t2) # 1991 1
ok <- ok & all(c(1991,1) == earliest.start(t1,t2))
ok

earliest.end(t1,t2)   # 1995 1
ok <- ok & all(c(1995,1) == earliest.end(t1,t2))
ok

latest.end(t1,t2)     # 1998 1
ok <- ok & all(c(1998,1) == latest.end(t1,t2))
ok

if (! ok) stop("some tests FAILED")

#  2000/03/21 14:42:28
#   A TS PADI server is necessary for these tests.
#   The next line is only necessary to remove this in an old version which set
#     home in frame 0. (I'll never do that again.)
#   if (is.S()) remove("DSE.HOME", where=0) 

   padi.function.tests.fame(verbose=F)
   simple.monitor.function.tests(verbose=F) 
   tagged.function.tests(verbose=F)
   combination.monitor.function.tests(verbose=F)
#  An RPC timeout sometimes happens in the fm tests above, but backward 
#      compatiability cannot support timeout.
   troll.function.tests(verbose=F)
# backward tests are not being kept up and may fail.
   fm.padi.function.tests(verbose=T) #backward compatabe FAME access using PADI

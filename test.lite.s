#   2000/03/21 14:42:25 

   system.info()
   version.dse()
   random.number.test() 
   tframe.function.tests(verbose=F) 
   dse1.function.tests(verbose=F)    
   example.verify.data(eg1.DSE.data.diff, fuzz.small=1e-12, verbose=F)  
# The following gave several..NOT CORRECT in R. It needs est.VARX.ar (not in earlier versions of R).
# example.tests(eg1.DSE.data.diff,fuzz.small=1e-12, verbose=T)
   dse2.function.tests(verbose=F, graphics=F)  
   dse3.function.tests(verbose=F, graphics=F) 
   dse4.function.tests(verbose=F, graphics=F) 
   dse2.graphics.tests(verbose=F)
   dse3.graphics.tests(verbose=F)
   dse4.graphics.tests(verbose=F)  #     test 3 needs stepwise
   guide.example.tests.part1(verbose=F, graphics=T)
#         gives  Warning cov. matrix is singular. Working on subspace

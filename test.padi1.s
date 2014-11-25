#   2000/03/21 14:42:28 
#######################################################################
#   Copyright 1993, 1994, 1995, 1996,  Bank of Canada.
#   Copyright 1997 (June), Paul Gilbert.
#   Copyright 1997 (Aug.), Bank of Canada.
#   Copyright 1998, Bank of Canada.

#   The user of this software has the right to use, reproduce and distribute it.
#   Bank of Canada makes no warranties with respect to the software or its 
#   fitness for any particular purpose. The software is distributed by the Bank
#   of Canada and by Paul Gilbert solely on an "as is" basis. By using the  
#   software, user agrees to accept the entire risk of using this software.

################################################################################

# These set of tests will not work unless the PADI interface is
#   also installed (available at http://www.bank-banque-canada.ca/pgilbert)

#    a TS PADI server is necessary for the following

   tfPADI.function.tests()
   padi.function.tests.simple(verbose=F)     # all ok
#   TSPADI.function.tests(verbose=T, ets=F)   # all ok, 
   TSPADI.function.tests(verbose=F, ets=T)   #  ets time out prob. sometimes

.First.lib <- function(library,section){
if(!require("dse1",  warn.conflicts=F)) 
    warning("This package requires the dse1 package.")
if(!require("padi",  warn.conflicts=F)) 
    warning("This package requires the padi package.")
invisible()}

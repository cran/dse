
###############################################

#   .tf and tftframe methods

# these provide .ts or tsp like methods for more general arrays with time as the
#   first dimension. It may be possible extend to more general structures.

# POSIX may be a better example

###############################################

#start.tf <- function(x) {start(tframe(x))}
#end.tf <- function(x) {end(tframe(x))}
#periods.tf <- function(x) {periods(tframe(x))}
#frequency.tf <- function(x) {frequency(tframe(x))}
#time.tf <- function(x) {time(tframe(x))}

start.tf <- function(x) {tframeStart(x)}
end.tf <- function(x) {tframeEnd(x)}
periods.tf <- function(x) {tframePeriods(x)}
frequency.tf <- function(x) {tframeFrequency(x)}
time.tf <- function(x) {tframeTime(x)}

tframeStart.tftframe <- function(x)
   {c(floor(x[1]), round(1 +(x[1]%%1)*x[3]))}

tframeEnd.tftframe <- function(x)
   {c(floor(x[2]), round(1 + (x[2]%%1)*x[3]))}

tframePeriods.tftframe <- function(x)  {1+round((x[2]-x[1])*x[3])}

tframeFrequency.tftframe <- function(x) {x[3]}

tframeTime.tftframe <- function(x) {x[1] + (seq(periods(x))-1)/x[3]}

tfwindow.tf <- function(x, start.=NULL, end.=NULL, warn=T, eps=.Options$ts.eps) 
  {# this needs work
   tfwindow(ts(x, start=start(x), frequency=frequency(x)),
             start.=start., end.=end., warn=warn)
  }


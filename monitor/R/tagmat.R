#   2000/04/18 11:15:54  
###########################################################################

# tagged data class  (matrix with a "tags" attribute)       <<<<<<<<<<<<

###########################################################################


##############################################################################

#  section containing documentation "stubs" (specific methods 
#  for generic functions) so that R CMD build does not complain.

##############################################################################



##############################################################################

#  end of section containing documentation "stubs" (specific methods 
#  for generic functions) so that R CMD build does not complain.

##############################################################################

tags <- function(x) {attr(x, "tags")}
 
"tags<-" <- function(x, value)   
  {if (is.null(value))
       {attr(x, "tags") <- NULL
        if (!is.null(class(x))) class(x) <- class(x)[ class(x) != "tagged"]
        return(x)
       }
   if (length(value) == 1) value <- array(value, dim(x))
    # drop any extra attributes
   attributes(value) <- list(dim=attributes(value)$dim)
   attr(x, "tags") <- value 
   classed(x, c("tagged", dseclass(x))) # constructor ("tags<-")
  }

tagged <- function(x, ...)UseMethod("tagged")

tagged.default <- function(x, tags) {tags(x) <- tags; x}


tagged.TSdata <- function(x, input.tags, output.tags)
  {if(0 != input.dimension(x))   input.data(tags(x)) <-  input.tags
   if(0 != output.dimension(x)) output.data(tags(x)) <- output.tags
   x
  }

select.series.tagged <- function(x, series=seq(ncol(x)))
     {names <- series.names(x)
      if (is.character(series)) series <- match(names,series, nomatch=0)
      tagged(select.series.default(x, series=series),
             select.series.default(tags(x), series=series))
     }

tbind.tagged <- function(mat1, mat2)
{# aline and bind ts matrices and tags
 if (is.tagged(mat1)) tag1 <- tags(mat1)
 else                 tag1 <- array("mat1", dim(mat1))
 if (is.tagged(mat2)) tag2 <- tags(mat2)
 else                 tag2 <- array("mat2", dim(mat2))
 tframe(tag1) <- tframe(mat1)
 tframe(tag2) <- tframe(mat2)
 cls <- dseclass(mat1)
 # this should use NextMethod
 dseclass(mat1) <- dseclass(mat1)[-1]  # otherwise tbind calls this tbind
 if (0 == length(dseclass(mat1))) dseclass(mat1) <- NULL
 dseclass(mat2) <- dseclass(mat2)[-1]  # otherwise tbind calls this tbind
 if (0 == length(class(mat2))) dseclass(mat2) <- NULL
 tagged(classed(tbind(mat1, mat2), cls),  tbind(tag1,tag2)) 
}

is.tagged <- function(obj)  {inherits(obj,"tagged")}

test.equal.tagged <- function(mat1, mat2)
{ test.equal.matrix(mat1,mat2) & 
  test.equal.matrix(attr(mat1,"tags"), attr(mat2, "tags"))
}



fprint <- function(matrix, super.title=NULL, sub.title=NULL, 
        digits=options()$digits, space=" ", file=NULL, append=F)
   {UseMethod("fprint")}

fprint.tagged <- function(matrix, super.title=NULL, sub.title=NULL, 
        digits=options()$digits, space=" ", file=NULL, append=F) 
 {# Formattted print of a matrix of class tagged.
  # Corresponding characters are printed after matrix numbers.
  # A character matrix (out) is returned invisibly.
  # If file is not NULL then elements of out are printed to lines of the file.
  tags <- attr(matrix, "tags")
  out <- NULL
  f <- frequency(matrix)
  s <- start(matrix)
  s <- s[1] + (s[2]-1)/f
  if (12 ==f) p <- c("Jan","Feb","Mar","Apr","May", "Jun","Jul","Aug", "Sep",
         "Oct","Nov","Dec")
  else if (4 == f) p <- c("Q1","Q2","Q3","Q4")
  else if (52 == f) p <- format(1:52)
  else p <-NULL
  pre.space <- paste(rep(" ",nchar(format(s))+nchar(p[1])),collapse="")
  if (!is.null(super.title))  out <- paste(pre.space, super.title, sep="")
  names <- format(dimnames(matrix)[[2]], digits=digits)
  if (!is.null(names))
    {ot <- pre.space
     for (i in seq(length(names)))
        ot <- paste(ot, space,names[i],sep="")
     out <- c(out, ot)
    }
  if (!is.null(sub.title)) out <- c(out,paste(pre.space, sub.title,sep=""))
  m <- format(signif(matrix[,], digits=digits))
  for (i in seq(nrow(m))) 
    {d <- (s+(i-1)/f) +.Options$ts.eps # +eps or trunc sometimes gets wrong year
     ot <- paste(trunc(d)," ", p[round(1+f*(d%%1))]," ", sep ="")
     for (j in seq(ncol(m))) 
       {ot <-paste(ot,space, m[i,j], sep="")
        if (!is.null(tags)) ot <- paste(ot,tags[i,j], sep="")
       }
      out <- c(out, ot)
    }
  if (!is.null(file)) write(out, file=file, append=append)
  invisible(out)
 }


splice.tagged <- function(mat1, mat2, tag1=tags(mat1), tag2=tags(mat2))
{# splice together 2 time series matrices as with splice.ts.
 # If data  is provided in both for a given period then mat1 takes priority.
 # The frequencies should be the same.
 # tag1 and tag2 are taken from mat1 and mat2 unless
 #   they are specified in the argument. If specified they
 #   should be single character strings or matrices of character 
 #   strings of same dimension as mat1 and mat2. This second is useful for multiple
 # applications of the function. The result is the
 # resulting spliced matrix of class "tagged"  
 # (suitable for use with fprint).
 # In the case tags are not available and are not specified 
 #   in the argument then they are set to "mat1" and "mat2".
 cls <- dseclass(mat1)
 if (is.null(tag1)) tag1 <- "mat1"
 if (is.null(tag2)) tag2 <- "mat2"
 if (length(tag1) == 1) tag1 <- array(tag1, dim(mat1))
 if (length(tag2) == 1) tag2 <- array(tag2, dim(mat2))
 if (is.null(mat1) & is.null(mat2)) return(NULL)
 if (is.null(mat2)) return(tagged(mat1, tag1))
 if (is.null(mat1)) return(tagged(mat2, tag2))
 freq <- frequency(mat1)
 if (freq != frequency(mat2)) stop("frequencies must be the same.\n")
 p <- dim(mat1)[2]
 if (p != dim(mat2)[2]) stop("number of series must be the same.\n")
 tframe(tag1) <- tframe(mat1)
 tframe(tag2) <- tframe(mat2)

 fr <- c(freq,1)
 st <- min(fr %*% start(mat1), fr %*% start(mat2))
 strt <- c(st %/% freq, st %% freq)
 en <- max(fr %*% end(mat1), fr%*% end(mat2))
 tf <- list(start=strt, frequency=freq)
 if (fr %*% start(mat1) > st) 
    {tag1 <-tframed(rbind(matrix("?", fr %*% start(mat1) -st, p), tag1),tf)
     mat1 <-tframed(rbind(matrix(NA,  fr %*% start(mat1) -st, p), mat1), tf)
    }
 if (fr %*%   end(mat1) < en) 
    {tag1 <-tframed(rbind(tag1, matrix("?", en - fr %*% end(mat1), p)), tf)
     mat1 <-tframed(rbind(mat1, matrix(NA,  en - fr %*% end(mat1), p)), tf)
    }
 if (fr %*% start(mat2) > st) 
    {tag2 <-tframed(rbind(matrix("?", fr %*% start(mat2) -st, p), tag2), tf)
     mat2 <-tframed(rbind(matrix(NA,  fr %*% start(mat2) -st, p), mat2), tf)
    }
 if (fr %*%   end(mat2) < en) 
    {tag2 <-tframed(rbind(tag2,matrix("?", en - fr %*% end(mat2), p)), tf)
     mat2 <-tframed(rbind(mat2, matrix(NA, en - fr %*% end(mat2), p)), tf)
    }
 na <- is.na(mat1)
#browser()
 mat1[na]  <- mat2[na]
 tag1[na] <- tag2[na]
 dimnames(mat1) <-list(round(time(mat1),digits=3),dimnames(mat1)[[2]])
 tags(mat1) <- tag1
 classed(mat1, cls )
}

trim.na.tagged <- function(mat, start.=T, end.=T)
{# trim NAs from the ends of a ts matrix of class "tagged".
 # (Observations for all series are dropped in a given period if any 
 #  one contains an NA in that period.)
 # if start.=F then beginning NAs are not trimmed.
 # If end.=F   then ending NAs are not trimmed.
 sample <- ! apply(is.na(mat),1, any)
 if (start.) s <-min(time(mat)[sample])
 else       s <-start(mat)
 if (end.)   e <-max(time(mat)[sample])
 else       e <-end(mat)
 tfwindow(mat,start=s, end=e, warn=F)
}

tfwindow.tagged <- function(x, start=NULL, end=NULL, warn=T)
{# window a ts matrix of class "tagged".
 # With the default warn=T warnings will be issued if no truncation takes
 #  place because start or end is outside the range of data.
 tags <- tags(x)
 dseclass(x) <- dseclass(x)[-1]
 if (0 == length(dseclass(x))) dseclass(x) <- NULL
 # The next line converts scalars tags to a matrix.
 if (length(tags) == 1) tags <- array(tags, dim(x))
 # The next lines converts missing tags to a matrix.
 if (length(tags) == 0)
   {tags <- array("", dim(x))
    if (warn) warning("missing tags converted to empty string.")
   }
 tframe(tags) <- tframe(x)
 # The following is complicated by the fact that some versions of window
 #    look for missing arguments.
 if (is.null(start))
   {x   <- tfwindow(x  , end=end, warn=warn)
    tags<- tfwindow(tags,end=end, warn=warn)
   }
 else if (is.null(end))
   {x   <- tfwindow(x  , start=start, warn=warn)
    tags<- tfwindow(tags,start=start, warn=warn)
   }
 else
   {x   <- tfwindow(x,   start=start, end=end, warn=warn)
    tags<- tfwindow(tags,start=start, end=end, warn=warn)
   }
 tagged(x, tags)
}




tagged.function.tests <- function(verbose=T, synopsis=T, fuzz.small=1e-10)
{# A short set of tests of the tagged class methods. 

  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (!is.TSdata(eg1.DSE.data.diff))
     stop("Test data not found. Testing stopped.")
  if (synopsis & !verbose) cat("All tagged class tests ...")
  if (verbose) cat("tagged class test 1 ... ")
#  z <- output.data(eg1.DSE.data.diff)
#  tags(z, "tags") <- array("a", dim(z))
#  dseclass(z) <- "tagged"
  z <- output.data(eg1.DSE.data.diff)
  z <- tagged(z, array("a", dim(z)))
  ok <- is.tagged(z)
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("tagged class test 2... ")
#  zz <- z
#  tags(zz) <- array("b", dim(z))
  zz <- tagged(z, array("b", dim(z)))
  ok <- test.equal(z,z) & (!test.equal(z,zz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 3... ")
  zz <- tfwindow(z, start=c(1989,1))
  tags(zz) <- array("b", dim(zz))
  zzz <- tbind(tfwindow(z, start=c(1989,1)),zz)
  ok <- (2*sum(tfwindow(output.data(eg1.DSE.data.diff),
           start=c(1989,1)))) ==  sum(zzz)
  ok <- ok & all("a" == tags(zzz)[,1:3]) &  all("b" == tags(zzz)[,4:6]) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 4... ")
  zzz <- splice(zz, tfwindow(z, end=c(1990,1)))
  ok <- test.equal.matrix(z,zzz) & (!test.equal(z,zzz))
  zzz <- splice(zz, tfwindow(output.data(eg1.DSE.data.diff),
                           end=c(1990,1)), tag2="x")
  ok <- ok & test.equal.matrix(z,zzz) & (!test.equal(z,zzz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All tagged class tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }
  invisible(all.ok)
}


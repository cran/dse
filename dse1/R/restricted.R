####################################################################
#     
#   Define a model with restrictions (or possibly extentions)
#    consider TSmodifiedModel
#     
####################################################################


TSrestrictedModel <- function(model,
                    coefficients=NULL,
                    restriction=NULL){
  # verify the structure
  if(!is.TSmodel(model))
       stop("The model should be a TSmodel.")
  if(inherits(model, "TSrestrictedModel"))
       stop("The model should not already be a TSrestrictedModel.")
  if(is.null(coefficients)) stop("A vector coefficients must be supplied.")
  if(is.null(restriction)) stop("A function restriction must be supplied.")
  # verify that the restriction works
  if(!is.TSmodel( restriction(model, coefficients)))
       stop("The restriction function is not returning a TSmodel.")
  if(inherits( restriction(model, coefficients), "TSrestrictedModel"))
       stop("The restriction function should not return a TSrestrictedModel.")
  classed(list(TSmodel=model,
               coefficients=coefficients, 
	       restriction=restriction),
       c("TSrestrictedModel","TSmodel")) # constructor
  }

coef.TSrestrictedModel <- function (object, ...){object$coefficients}

setArrays.TSrestrictedModel <- function(model, coefficients=coef(model))
  {#coefficients is coefficients (restriction coefs appended 
   #  with coef for model$TSmodel )
   # This does not really set arrays, but relies on model$restricton to be used
   # in any evaluation of the model, and it must set arrays using coefficients
   model$coefficients <- coefficients
   model
   }

setTSmodelParameters.TSrestrictedModel <- function(model,
         constants=TSmodel(model$TSmodel)$constants) { 
  # not sure what is really needed her
  #setTSmodelParameters(TSmodel(model), constants=constants)
  model
  }

print.TSrestrictedModel <- function (x, ...){
   print(x$restriction(x$TSmodel, x$coefficients), ...)
   }

l.TSrestrictedModel <- function (obj1, obj2, result=NULL, ...){
   r <- l(obj1$restriction(obj1$TSmodel, obj1$coefficients), obj2,
          result=result, ...)
   if ( (!is.null(result)) && (result =="like")) return(r) # neg.log.like. from residualStats
   r$model <- TSrestrictedModel(r$model, 
                      coefficients=obj1$coefficients,
                      restriction=obj1$restriction)

   if ( is.null(result)) r else r[[result]] 
   }


smoother.TSrestrictedModel <- function(model, data, compiled=.DSEflags()$COMPILED){
    m <- model$restriction(model$TSmodel, model$coefficients)
    r <- smoother(m, data, compiled=compiled)
    r$model <- model$TSmodel
    r
    }

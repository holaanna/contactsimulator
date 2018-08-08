## basic methods
setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("window",function(x,...)standardGeneric("window"))

## constituent components of a 'contact' object
setGeneric("r.model",function(object,...)standardGeneric("r.model"))
setGeneric("d.model",function(object,...)standardGeneric("d.model"))
setGeneric("rmodel",function(object,...)standardGeneric("rmodel"))
setGeneric("kernelmodel",function(object,...)standardGeneric("kernelmodel"))
setGeneric("d.prior",function(object,...)standardGeneric("d.prior"))
setGeneric("r.prior",function(object,...)standardGeneric("r.prior"))
setGeneric("init.state",function(object,...)standardGeneric("init.state"))

## internals
# setGeneric("contact.fun",function(f,...)standardGeneric("contact.fun"))
## functions to extract or call the components of a "pomp" object
setGeneric("obs",function(object,...)standardGeneric("obs"))
setGeneric("time",function(x,...)standardGeneric("time"))
setGeneric("time<-",function(object,...,value)standardGeneric("time<-"))
setGeneric("coef",function(object,...)standardGeneric("coef"))
setGeneric("coef<-",function(object,...,value)standardGeneric("coef<-"))
setGeneric("states",function(object,...)standardGeneric("states"))
setGeneric("timezero",function(object,...)standardGeneric("timezero"))
setGeneric("timezero<-",function(object,...,value)standardGeneric("timezero<-"))
#setGeneric("partrans",function(object,params,dir,...)standardGeneric("partrans"))
#setGeneric("logLik",function(object,...)standardGeneric("logLik"))


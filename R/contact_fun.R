#============================================================================
# A class for functions that may be defined in R, using pre-written native routines of Cpp functions
#============================================================================
#' A class for functions that may be defined in R, using pre-written native routines of Cpp functions
#'
#' Definition and methods of the \sQuote{contact_fun} class.
#'
#' The \sQuote{contact_fun} class implements a common interface for user-defined procedures that can be defined in terms of R code or by compiled native routines.
#'
#' @name contact_fun
#' @rdname contact_fun
#' @include contactsimulator-package.R cppsnippet.R rsnippet.R generics.R
#' @docType methods
#' @keywords internal
#' @concept extending the contact package
#' @concept low-level interface
#'
#' @param f A function or the name of a native routine.
#' @param PACKAGE optional; the name of the dynamically-loadable library in
#' which the native function \code{f} can be found.
#' @param proto optional string; a prototype against which \code{f} will be
#' checked.
#' @param object,x the \sQuote{contact_fun} object.
#' @author Hola K. Adrakey
#' @seealso \code{\link{contact}}
NULL

contactfunmode <- list(undef=-1L, Rfun=0L, native=1L, regNative=2L)


setClass(
  "contact_fun",
  slots = c(
    R.fun = "function",
    native.fun = "character",
    PACKAGE = "character",
    mode = "integer",
    address = "externalptr",
    obsnames = 'character',
    statenames = 'character',
    paramnames = 'character',
    covarnames = 'character',
    purpose = "character"
  ),
  prototype = prototype(
    R.fun=function (...) {
      stop("contact_fun","unreachable error: please report this bug!")
    },
    native.fun=character(0),
    PACKAGE = character(0),
    mode =  contactfunmode$undef,
    obsnames = character(0),
    statenames = character(0),
    paramnames = character(0),
    covarnames = character(0),
    purpose = "a needed function"
  )

)

setGeneric("contact_fun", function(f, ...) standardGeneric("contact_fun"))


setMethod(
  "contact_fun",
  signature = signature(f="missing"),
  definition = function(f, slotname = NULL, ...){
    new("contact_fun", purpose=as.character(slotname))
  }
  )


setMethod(
  "contact_fun",
  signature = signature(f="NULL"),
  definition = function(f, slotname = NULL, ...){
    new("contact_fun", purpose=as.character(slotname))
  }
)

setMethod(
  "contact_fun",
  signature = signature(f="ANY"),
  definition=function (f, slotname = NULL, ...) {
    stop("bad option for ",sQuote(slotname)," argument",call.=FALSE)
  }
)


setMethod(
  "contact_fun",
  signature=signature(f="function"),
  definition=function (f, proto = NULL, slotname = NULL, ...) {
    if (!is.null(proto)) {
      prototype <- as.character(proto)
      fname <- prototype[1]
      args <- prototype[-1]
      if (is.function(f)&&(!all(args%in%names(formals(f)))))
        stop(slotname,
          sQuote(fname)," must be a function of prototype ",
          sQuote(deparse(proto)))
    }
    new("contact_fun",R.fun=f,mode=pompfunmode$Rfun,purpose=as.character(slotname))
  }
)


setMethod(
  "contact_fun",
  signature=signature(f="Rsnippet"),
  definition=function (f, PACKAGE, slotname = NULL, libname = NULL,
                       obsnames = character(0), statenames = character(0),
                       paramnames = character(0), covarnames = character(0), ...) {
    if (is.null(slotname))
      stop("in ",sQuote("contact_funn"),": unspecified ",
           sQuote("slotname"),call.=FALSE)
    if (is.null(libname))
      stop("in ",sQuote("contact_fun"),": unspecified ",
           sQuote("libname"),call.=FALSE)
    slotname <- as.character(slotname)
    libname <- as.character(libname)
    new(
      "contact_fun",
      native.fun=render(slotname,name=libname),
      PACKAGE=as.character(libname),
      mode=contactfunmode$Rfun,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=slotname
    )
  }
)


setMethod(
  "contact_fun",
  signature=signature(f="Cppsnippet"),
  definition=function (f, PACKAGE, slotname = NULL, libname = NULL,
                       obsnames = character(0), statenames = character(0),
                       paramnames = character(0), covarnames = character(0), ...) {
    if (is.null(slotname))
      stop("in ",sQuote("contact_fun"),": unspecified ",
           sQuote("slotname"),call.=FALSE)
    if (is.null(libname))
      stop("in ",sQuote("contact_fun"),": unspecified ",
           sQuote("libname"),call.=FALSE)
    slotname <- as.character(slotname)
    libname <- as.character(libname)
    new(
      "contact_fun",
      native.fun=render(slotname,name=libname),
      PACKAGE=as.character(libname),
      mode=contactfunmode$regNative,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=slotname
    )
  }
)



setMethod(
  "contact_fun",
  signature = signature(f="character"),
  definition = function(f, PACKAGE = NULL,
                        slotname = NULL, ...){
    new(
      "contact_fun",
      native.fun = f,
      PACKAGE = as.character(PACKAGE),
      mode = contactfunmode$native,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=as.character(slotname)
    )
  }
)

setMethod(
  "contact_fun",
  signature=signature(f="contact_fun"),
  definition=function (f, ...) f
)

setMethod(
  'show',
  signature=signature('contact_fun'),
  definition=function (object) {
    mode <- object@mode
    if (mode==contactfunmode$Rfun) { # R function
      show(object@R.fun)
    } else if (mode==contactfunmode$native) { # user supplied native code
      cat("native function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
      cat ("\n")
    } else if (mode==contactfunmode$regNative) { # built from C snippets
      cat("native function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", defined by a Cppsnippet",sep="")
      cat ("\n")
    } else {
      cat("not specified\n")
    }
  }
)

setMethod(
  'print',
  'contact_fun',
  function (x, ...) show(x)
)

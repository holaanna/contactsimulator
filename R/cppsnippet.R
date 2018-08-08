## a class to hold snippets of C code

setClass(
  "Cppsnippet",
  slots=c(
    text="character"
  ),
  prototype=prototype(
    text=character(0)
  )
)

Cppsnippet <- function (text) {
  new(
    "Cppsnippet",
    text=as.character(text)
  )
}

setAs(
  from="Cppsnippet",
  to="character",
  def = function (from) {
    from@text
  }
)

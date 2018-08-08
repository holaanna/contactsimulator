## a class to hold snippets of C code

setClass(
  "Rsnippet",
  slots=c(
    text="character"
  ),
  prototype=prototype(
    text=character(0)
  )
)

Rsnippet <- function (text) {
  new(
    "Rsnippet",
    text=as.character(text)
  )
}

setAs(
  from="Rsnippet",
  to="character",
  def = function (from) {
    from@text
  }
)

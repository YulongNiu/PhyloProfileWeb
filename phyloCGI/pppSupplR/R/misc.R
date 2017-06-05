##' Miscellaneous - Retreive grepexpr content
##'
##' Retreive the content from the output of grepexpr() function.
##' @title Retreive grepexpr content
##' @param s Input string.
##' @param g A list that output from grepexpr()
##' @return A string
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##' 
getcontent <- function(s, g) {
  substring(s, g, g+attr(g,'match.length')-1)
}

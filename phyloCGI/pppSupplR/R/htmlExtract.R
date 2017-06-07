##' Extract html table content from a vector
##'
##' Extract the <tr>...</tr> from a vector which is yield from the "hwriter" pacakge. 
##' @title Extract html table 
##' @param tableVec A vector containg a entril html table.
##' @return A vector only containing <tr>...</tr>
##' @examples
##' htmlVec <- c('<table border=\"1\">\n<tr>\n<td>1</td><td>3</td>
##' </tr>\n<tr>\n<td>2</td><td>4</td></tr>\n</table>\n')
##' tableVec <- htmlExtractTable(htmlVec)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
htmlExtractTable <- function(tableVec) {

  trReg <- gregexpr('<tr>.*</tr>', tableVec)
  trVal <- getcontent(tableVec, trReg[[1]])

  return(trVal)
}

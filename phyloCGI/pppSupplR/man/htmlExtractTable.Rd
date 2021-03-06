% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htmlExtract.R
\name{htmlExtractTable}
\alias{htmlExtractTable}
\title{Extract html table}
\usage{
htmlExtractTable(tableVec)
}
\arguments{
\item{tableVec}{A vector containg a entril html table.}
}
\value{
A vector only containing <tr>...</tr>
}
\description{
Extract html table content from a vector
}
\details{
Extract the <tr>...</tr> from a vector which is yield from the "hwriter" pacakge.
}
\examples{
htmlVec <- c('<table border=\\"1\\">\\n<tr>\\n<td>1</td><td>3</td>
</tr>\\n<tr>\\n<td>2</td><td>4</td></tr>\\n</table>\\n')
tableVec <- htmlExtractTable(htmlVec)
}
\author{
Yulong Niu \email{niuylscu@gmail.com}
}

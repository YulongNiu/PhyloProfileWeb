GetLinkages <- function(geneIDs, linkData) {

  ## USE: get the linkages from 'linkData'
  ## INPUT: 'geneIDs' is the gene ids. 'linkData' is a matrix, of which the first and third columns is gene Ids.
  ## OUTPU: linkage matrix

  fromIdx <- linkData[, 1] %in% geneIDs
  toIdx <- linkData[, 3] %in% geneIDs

  linkIdx <- fromIdx | toIdx
  candLinksMat <- linkData[linkIdx, , drop = FALSE]

  return(candLinksMat)
  
}

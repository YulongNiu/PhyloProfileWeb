##' Generate phylogenetic data for Circos plot
##'
##' This function is used to generate Circos links, Circos lables, and Circos heatmap.
##' @title Generate Circos phylogenetic files
##' @param ft A from-to matrix, here we used the undirected network. The 1st and 2nd columns are gene ID, and the 3rd clumn is the weight ranging from 0 to 1 like correlation coefficient or similarity value.
##' @param locaAnno The gene location annotation information. The 1st column is the gene ID which should be the same format in "ft". The 2nd column is the chromosome name, like "1", "2", and "MT". The 3rd and 4th columns are the start and end position, respectively.
##' @param nodeName Setting node names desided to show. The default value is "all" that showing all the input nodes. Users can set this parameters to a subset of the nodes.
##' @param thick The thickness
##' @param showEdge Whether or not to show the weighted thickness (wThickness)of linkages. The default value is "TRUE" that means the wThickness = thick * weight. Otherwise, the thickness of all the linkages is set the same as the input "thick".
##' @return A list of Circos links, Circos lables, and Circos heatmap.
##' @examples
##' data(atpft)
##' data(geneAnno)
##' # preprocess
##' geneAnno <- geneAnno[, -2]
##' interLogic <- (atpft[, 1] %in% geneAnno[, 1]) & (atpft[, 2] %in% geneAnno[, 1])
##' atpft <- atpft[interLogic, ]
##' # show all the linkages with weighted thickness
##' fatpCircos <- ft2circos(ft = atpft, locaAnno = geneAnno)
##' # only show the links with "ATP5A1" without weighed thickness
##' fatpCircos <- ft2circos(ft = atpft, locaAnno = geneAnno, nodeName = 'ATP5A1', showEdge = FALSE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
ft2circos <- function(ft, locaAnno, nodeName = 'all', thick = 5, showEdge = TRUE){

  if (!identical(nodeName, 'all')){
    ftIn <- ft[ft[, 1] %in% nodeName, , drop = FALSE]
    ftOut <- ft[ft[, 2] %in% nodeName, , drop = FALSE]
    ftOut[, 1:2] <- ftOut[, 2:1]
    ft <- rbind(ftIn, ftOut)
  } else {}

  locaAnno[, 2] <- paste('hs', locaAnno[, 2], sep = '')

  # get nodes genome location information
  weight <- as.numeric(as.character(ft[, 3]))
  nodeIn <- locaAnno[match(ft[, 1], locaAnno[, 1]), 1:4]
  nodeOut <- locaAnno[match(ft[, 2], locaAnno[, 1]), 1:4]

  # get gene lables data
  circosLabel <- rbind(nodeIn[, c(2:4, 1)], nodeOut[, c(2:4, 1)])
  circosLabel <- circosLabel[!duplicated(circosLabel, MARGIN = 1), ]
  rownames(circosLabel) <- NULL

  # get link data
  circosLink <- cbind(nodeIn[, 2:4], nodeOut[, 2:4])
  rownames(circosLink) <- NULL

  # show edge
  # in case of no nodes and edges
  if (dim(circosLink)[1] != 0){
    if (showEdge){
      weight <- weight * thick
      weight <- paste('thickness=', weight, 'p', sep = '')
      circosLink <- cbind(circosLink, weight)
    } else {
      circosLink <- cbind(circosLink, paste('thickness=', thick, 'p', sep = ''))
    }
  } else {}

  # show heatmap
  if (dim(circosLink)[1] != 0){
    heatVal <- cbind(circosLink[, 4:6], weight)
    ## heatVal <- apply(heatVal, 1:2, as.character)
    ## heatVal <- rbind(heatVal, c(as.character(circosLink[1, 1:3]), 1))
  } else {
    heatVal <- circosLink[, 1:4, drop = FALSE]
  }

  return(list(circosLink = circosLink, circosLabel = circosLabel, circosHeatmap = heatVal))

}



##' Generate gene frequency from phylogenetic profiles.
##'
##' Calculate gene frequence of three domains of life or more details of Eukaryota.
##' @title Generate gene frequency
##' @param KEGGPhylo  A character matrix, the 1st column is the KEGGID, and the 2nd column is species phylogeny for example "Eukaryotes;Animals;Vertebrates;Mammals".
##' @param selectPhyloData  A numeric matrix. The row names are the genes and and the column names are the species.
##' @param splitEu Whether to show the eukaryotic split, the default value is "FALSE".
##' @return A fequency matrix.
##' @examples
##' data(phyloSpe)
##' data(wholeProfile)
##' freqMat <- SpeFreq(KEGGPhylo = phyloSpe, selectPhyloData = wholeProfile, splitEu = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
SpeFreq <- function(KEGGPhylo, selectPhyloData, splitEu = FALSE) {

  phyloCode <- as.character(KEGGPhylo[, 2])
  phyloCode <- sapply(strsplit(phyloCode, split = ';', fixed = TRUE), function(x) x[1:2])
  phyloCode <- t(phyloCode)
  phyloCode <- cbind(as.character(KEGGPhylo[, 1]), phyloCode)

  phyloList <- list()
  phyloList$arcSpe <- phyloCode[phyloCode[, 3] %in% 'Archaea', 1]
  phyloList$bacSpe <- phyloCode[phyloCode[, 3] %in% 'Bacteria', 1]
  if (!splitEu) {
    phyloList$euSpe <- phyloCode[!(phyloCode[, 3] %in% c('Archaea', 'Bacteria')), 1]
  } else {
    phyloList$Animals <- phyloCode[phyloCode[, 3] %in% 'Animals', 1]
    phyloList$Plants <- phyloCode[phyloCode[, 3] %in% 'Plants', 1]
    phyloList$Fungi <- phyloCode[phyloCode[, 3] %in% 'Fungi', 1]
    phyloList$Protists <- phyloCode[phyloCode[, 3] %in% 'Protists', 1]
  }

  phyloFreqList <- lapply(phyloList, function(x) {
    phyloEachData <- selectPhyloData[, colnames(selectPhyloData) %in% x]
    phyloEachFreq <- apply(phyloEachData, 1, function(x) {sum(x)/length(x)})
    return(phyloEachFreq)
  })

  phyloFreqMat <- do.call(cbind, phyloFreqList)

  # asign names
  if (!splitEu) {
    colnames(phyloFreqMat) <- c('freqArc', 'freqBac', 'freqEu')
  } else {
    colnames(phyloFreqMat) <- c('freqArc', 'freqBac', 'freqAni', 'freqPlant', 'freqFun', 'freqProt')
  }

  return(phyloFreqMat)
}




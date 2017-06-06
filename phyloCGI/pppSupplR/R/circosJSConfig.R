##' Generate gene frequency for circosJS from phylogenetic profiles.
##'
##' Calculate gene frequence of three domains of life or more details of Eukaryota.
##'
##' @title Generate gene frequency for circosJS
##' @param allAnno The annotation matrix. The 1st is the KEGG gene IDs, which has the same format with "inputft". The 2nd is the symbol name shown in the Circos plot. The 3rd column is the chromosome name. The 4th and 5th columns are the start and end positions. The "allAnno" may contains genes that are absent in the "inputft".
##' @param allSpePhylo The species and phylogeny matrix.
##' @param allPhyloData The phylogenetic profiles.
##' @return A fequency matrix.
##' @examples
##' data(phyloSpe)
##' data(wholeProfile)
##' data(geneAnno)
##' freqMat <- SpeFreqJS(phyloSpe, wholeProfile, geneAnno, splitEu = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
SpeFreqJS <- function(allSpePhylo, allPhyloData, allAnno, splitEu = FALSE) {

  ## step1 select profiles of genes with annotation
  allPhyloData <- allPhyloData[rownames(allPhyloData) %in% allAnno[, 1], , drop = FALSE]

  ## step2 split phylogeny
  phyloCode <- as.character(allSpePhylo[, 2])
  phyloCode <- sapply(strsplit(phyloCode, split = ';', fixed = TRUE), function(x) x[1:2])
  phyloCode <- t(phyloCode)
  phyloCode <- cbind(as.character(allSpePhylo[, 1]), phyloCode)

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

  ## step3 calculate freq
  phyloFreqList <- lapply(phyloList, function(x) {
    phyloEachData <- allPhyloData[, colnames(allPhyloData) %in% x]
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

  ## step4 merge anno
  anno <- allAnno[match(rownames(phyloFreqMat), allAnno[, 1]), -1]
  anno[, 2] <- paste0('chr', anno[, 2])
  phyloFreqMat <- cbind(phyloFreqMat, anno)

  return(phyloFreqMat)
}

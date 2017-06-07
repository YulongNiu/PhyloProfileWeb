##' Generate gene frequency for circosJS from phylogenetic profiles.
##'
##' Calculate gene frequence of three domains of life or more details of Eukaryota.
##'
##' @title Generate gene frequency for circosJS
##' @param allSpePhylo The species and phylogeny matrix.
##' @param allPhyloData The phylogenetic profiles.
##' @param splitEu Whether to show the eukaryotic split, the default value is "FALSE".
##' @inheritParams CheckLinkCol
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
  colnames(anno) <- c('label', 'chromosome', 'start', 'end')
  phyloFreqMat <- cbind(phyloFreqMat, anno)

  return(phyloFreqMat)
}


##' Generate linkages for circosJS plot
##'
##' This function is used to generate linkages for circosJS plot
##' @title Linkages for circosJS plot
##' @param ft A from-to matrix, here we used the undirected network. The 1st and 2nd columns are gene ID, and other columns are similarity or distances (Jaccard, Cor et. al.).
##' @inheritParams SpeFreqJS
##' @return A linkage matrix.
##' @examples
##' data(atpft)
##' atpft <- atpft[, c(1, 3, 5:6)]
##' data(geneAnno)
##' LinkJS(atpft, geneAnno)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
LinkJS <- function(ft, allAnno) {
  ## step 1 select genes with anno
  isAnno <- (ft[, 1] %in% allAnno[, 1]) & (ft[, 2] %in% allAnno[, 1])
  ft <- ft[isAnno, ]
  colnames(ft)[3:4] <- c('Jaccard', 'Cor')

  ## step 2 anno
  fromAnno <- allAnno[match(ft[, 1], allAnno[, 1]), c(3:5, 2)]
  fromAnno[, 1] <- paste0('chr', fromAnno[, 1])
  colnames(fromAnno) <- c('source_id', 'source_start', 'source_end', 'source_label')
  rownames(fromAnno) <- NULL
  toAnno <- allAnno[match(ft[, 2], allAnno[, 1]), c(3:5, 2)]
  toAnno[, 1] <- paste0('chr', toAnno[, 1])
  colnames(toAnno) <- c('target_id', 'target_start', 'target_end', 'target_label')
  rownames(toAnno) <- NULL

  linkMat <- cbind(fromAnno,
                   toAnno,
                   ft[, 3:4])

  return(linkMat)

}



##' Automatically write circosJS data files for plotting linkages
##'
##' Generate two types circosJS data files:
##' 1. linkages for each input "geneVec";
##' 2. frequency of presence in three domains;
##
##' @title Write circosJS data files
##' @param geneVecList The output object from the "CheckLinkCol" function.
##' @param savePath The path used to store circosJS files.
##' @inheritParams SpeFreqJS
##' @inheritParams LinkJS
##' @return NULL
##' @examples
##' \dontrun{
##' data(atpft)
##' atpft <- atpft[, c(1, 3, 5:6)]
##' data(geneAnno)
##' data(phyloSpe)
##' data(wholeProfile)
##' g1 <- c('hsa:111111', 'hsa:498', 'hsa:506', 'hsa:509', 'hsa:516', 'hsa:517')
##' g1Col <- c('#36A03F', '#49AE4F', '#1C386A', '#A03531', NA, NA)
##' checkList <- CheckLinkCol(g1, g1Col, geneAnno)
##' tmp1 <- writeCircosJS(checkList, atpft, geneAnno, phyloSpe, wholeProfile, savePath = 'tmp1/')}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom utils write.csv
##' @export
##'
writeCircosJS <- function(geneVecList, ft, allAnno, allSpePhylo, allPhyloData, savePath){

  geneVec <- geneVecList[[1]]
  geneCol <- geneVecList[[2]]
  geneSym <- geneVecList[[3]]

  ## step 1 remove linkages do not contain geneVec
  has1 <- ft[, 1] %in% geneVec
  ftleft <- ft[!has1, , drop = FALSE]
  ftleft <- ftleft[ftleft[, 2] %in% geneVec, c(2:1, 3:4), drop = FALSE]
  colnames(ftleft)[1:2] <- colnames(ftleft)[2:1]
  ft <- rbind(ft[has1, , drop = FALSE],
              ftleft)

  ## step2 linkages
  links <- LinkJS(ft, allAnno)
  ## with colors
  links <- cbind(links,
                 linkage_color = geneCol[match(links[, 'source_label'], geneSym)])
  write.csv(links,
            paste0(savePath, 'linkage.csv'),
            row.names = FALSE,
            quote = FALSE)

  ## step3 frequence
  genes <- unique(c(ft[,1 ], ft[, 2]))
  freqGene <- SpeFreqJS(allSpePhylo,
                        allPhyloData[rownames(allPhyloData) %in% genes, , drop = FALSE],
                        allAnno,
                        splitEu = FALSE)
  for (i in 1:3) {
    freqMat <- freqGene[, c(i, 4:7)]
    colnames(freqMat)[1] <- 'value'
    write.csv(freqMat,
              paste0(savePath, colnames(freqGene)[i], '.csv'),
              row.names = FALSE,
              quote = FALSE)
  }

  return(NULL)
}


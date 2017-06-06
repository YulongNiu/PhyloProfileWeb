##' Automatically write circos data files for plotting linkages
##'
##' Generate four types circos data files:
##' 1. linkages for each input "geneVec";
##' 2. "heatmap"(correlation value) for each linkage;
##' 3. frequency of presence in three domains;
##' 4. names and locations of genes
##
##' @title Write circos data files
##' @param geneVec A vector containing condidate genes used for Circos plot. The length of one is also allowed.
##' @param inputft The "from-to" matrix. The 1st and 2nd columns are nodes names, and the 3rd column is the similarity. The "inputft" may contain genes not presenting the "geneVec".
##' @param allAnno The annotation matrix. The 1st is the KEGG gene IDs, which has the same format with "inputft". The 2nd is the symbol name shown in the Circos plot. The 3rd column is the chromosome name. The 4th and 5th columns are the start and end positions. The "allAnno" may contains genes that are absent in the "inputft".
##' @param allSpePhylo The species and phylogeny matrix.
##' @param allPhyloData The phylogenetic profiles.
##' @param savePath The path used to store Circos files.
##' @return NULL
##' @examples
##' \dontrun{
##' data(atpft)
##' data(geneAnno)
##' data(phyloSpe)
##' data(wholeProfile)
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' f0genes <- c('hsa:516', 'hsa:517', 'hsa:518', 'hsa:515', 'hsa:521',
##' 'hsa:522', 'hsa:9551', 'hsa:10476', 'hsa:10632', 'hsa:4508', 'hsa:4509')
##' fgenes <- c(f1genes, f0genes)
##' tmp1 <- writeCircos(f1genes, atpft, geneAnno, phyloSpe, wholeProfile, savePath = 'tmp1/')}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom utils write.table
##' @export
##' 
writeCircos <- function(geneVec, inputft, allAnno, allSpePhylo, allPhyloData, savePath){

  ## select ft matrix containing the 'geneVec'
  ftLogic <- (inputft[, 1] %in% geneVec) | (inputft[, 2] %in% geneVec)
  inputft <- inputft[ftLogic, ]

  ## select anno matrix
  annoLogic <- (inputft[, 1] %in% allAnno[, 1]) & (inputft[, 2] %in% allAnno[, 1])
  inputft <- inputft[annoLogic, ]

  ## names and locations of genes
  labelGene <- ft2circos(inputft, allAnno[, -2], geneVec)[[2]]
  labelGeneVec <- labelGene[, 4]
  labelGene[, 4] <- allAnno[match(labelGene[, 4], allAnno[, 1]), 2]
  write.table(labelGene, paste0(savePath, 'labelGene.txt'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

  ## linkages and heatmaps
  for (i in 1:length(geneVec)){
    corPhyloCir <- ft2circos(inputft, allAnno[, -2], geneVec[i], showEdge = FALSE, thick = 3)
    write.table(corPhyloCir[[1]], paste0(savePath, geneVec[i], '.txt'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(corPhyloCir[[3]], paste0(savePath, geneVec[i], 'heatmap', '.txt'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  ## frequency
  freqGene <- SpeFreq(allSpePhylo, allPhyloData)
  freqGene <- freqGene[match(labelGeneVec, rownames(freqGene)), ]

  for (i in 1:ncol(freqGene)) {
    freqMat <- cbind(labelGene[, 1:3], freqGene[, i])
    write.table(freqMat, paste0(savePath, colnames(freqGene)[i], '.txt'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  return(NULL)
}

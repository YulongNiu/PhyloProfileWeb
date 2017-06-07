##' Retrieve top linkages
##'
##' Give users freedom to choose threshold of top linkages.
##' @title Retrieve top linkages data
##' @param geneIDs The vector of input geneIDs.
##' @param linkData A list, first element is the n x m matrix, n is the top jaccard value (decreasing order), m is the whole genome genes. Second and third element is n x m matrix with Jaccard and Correlation values.
##' @param annoVec A vector contains m anntoations of m genes.
##' @param threshold A threshold should 0 < threshold <= n
##' @return a formated linkages
##' @examples
##' lkName <- matrix(ncol = 10, nrow = 9)
##' colnames(lkName) <- paste0('g', 1:10)
##' for(i in 1:10) {
##'   rc <- sample(1:10, 10)
##'   rc <- rc[rc != i]
##'   rcn <- paste0('g', rc)
##'   lkName[, i] <- rcn
##' }
##' lkJac <- matrix(round(runif(9 * 10), 2), ncol = 10)
##' colnames(lkJac) <- colnames(lkName)
##' lkCor <- matrix(round(runif(9 * 10),2),  ncol = 10)
##' colnames(lkCor) <- colnames(lkName)
##' lkData <- list(lkName, lkJac, lkCor)
##' gnIDs <- paste0('g', c(3, 5, 7))
##' aoVec <- paste0('gene', 1:10)
##' names(aoVec) <- colnames(lkName)
##' canLink <- GetTopLink(gnIDs, lkData, aoVec, 4)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
GetTopLink <- function(geneIDs, linkData, annoVec, threshold = 400) {
  
  ## select data based on threshold
  lkThresList <- lapply(linkData, function(x) {
    return(x[1:threshold, , drop = FALSE])
  })

  ## whole gene list
  lkNameMat <- lkThresList[[1]]
  geneNames <- colnames(lkNameMat)
  geneNum <- length(geneNames)

  ## filter matrix
  filterMat <- matrix(ncol = geneNum,
                      nrow = threshold)
  
  for (i in 1:geneNum) {
    if (geneNames[i] %in% geneIDs) {
      filterMat[, i] <- TRUE
    } else {
      filterMat[, i] <- lkNameMat[, i] %in% geneIDs
    }
  }

  ftMat <- cbind(rep(geneNames, apply(filterMat, 2, sum)),
                 lkNameMat[filterMat])


  jacCorMat <- cbind(c(lkThresList[[2]][filterMat]),
                     c(lkThresList[[3]][filterMat]))

  ## delet repeat
  bigIdx <- ftMat[, 1] <= ftMat[, 2]
  ftMat[bigIdx, ] <- ftMat[bigIdx, 2:1]
  dupIdx <- duplicated(ftMat)

  ftMat <- ftMat[!dupIdx, ]
  jacCorMat <- jacCorMat[!dupIdx, ]

  ## annotation
  fromAnno <- annoVec[match(ftMat[, 1], names(annoVec))]
  toAnno <- annoVec[match(ftMat[, 2], names(annoVec))]

  ## summary
  candLinksMat <- cbind(ftMat[, 1], fromAnno, ftMat[, 2], toAnno, jacCorMat)
  colnames(candLinksMat) <- c('From', 'FromAnno', 'To', 'ToAnno', 'Jaccard', 'Cor')
  rownames(candLinksMat) <- NULL

  ## sort and move input genes left
  adIdx <- !(candLinksMat[, 1] %in% geneIDs)
  candLinksMat[adIdx, 1:4] <- candLinksMat[adIdx, c(3:4, 1:2)]

  ## order by Jaccard and InputName
  candLinksMat <- candLinksMat[order(as.numeric(candLinksMat[, 5]), decreasing = TRUE), ]
  candLinksMat <- candLinksMat[order(candLinksMat[, 1]), ]

  return(candLinksMat)
  
}

## ##' Retrieve linkages 
## ##'
## ##' Retrieve the linkages containing the input gene list. This function will be replaced by a SQL query.
## ##' @title Retrieve linkages data
## ##' @param geneIDs The vector of geneIDs.
## ##' @param linkData A matrix, of which the first and third columns is gene Ids.
## ##' @return A linkage matrix.
## ##' @examples
## ##' genes <- c('a', 'c')
## ##' linkMat <- matrix(c(letters[1:5], 1:5, letters[5:1], 5:1),
## ##' ncol = 4,
## ##' nrow = 5,
## ##' dimnames = list(paste0('link', 1:5), c('From', 'FromLink', 'To', 'ToLink')))
## ##' geneLinkMat <- GetLinkages(genes, linkMat)
## ##' @author Yulong Niu \email{niuylscu@@gmail.com}
## ##' @export
## ##' 
## GetLinkages <- function(geneIDs, linkData) {

##   fromIdx <- linkData[, 1] %in% geneIDs
##   toIdx <- linkData[, 3] %in% geneIDs

##   linkIdx <- fromIdx | toIdx
##   candLinksMat <- linkData[linkIdx, , drop = FALSE]

##   ## sort and move input genes left
##   adIdx <- !(candLinksMat[, 1] %in% geneIDs)
##   candLinksMat[adIdx, 1:4] <- candLinksMat[adIdx, c(3:4, 1:2)]

##   candLinksMat <- candLinksMat[order(candLinksMat[, 1]), ]

##   return(candLinksMat)
  
## }

##' Retrieve profiles
##'
##' Retrieve the phylogenetic profiles containing the input gene list. This function will be replaced by a SQL query.
##' @title Retrieve profiles
##' @param profileData A numeric matrix. The row names are the genes and and the column names are the species.
##' @param geneIDs The vector of geneIDs.
##' @return The selected profiles.
##' @examples
##' genes <- c('a', 'c')
##' profileMat <- matrix(sample(0:1, size = 20, replace = TRUE),
##' ncol = 4,
##' nrow = 5,
##' dimnames = list(letters[1:5], paste0('spe', 1:4)))
##' geneProfileMat <- GetProfile(genes, profileMat)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
GetProfile <- function(geneIDs, profileData) {

  candProfile <- profileData[rownames(profileData) %in% geneIDs, , drop = FALSE]

  return(candProfile)
}

##' Check the input genes whether or not have annotation information and whether or not have colors
##'
##' Check the input genes with annotation or having colors. if not, an error message with length 1 will be returned.
##' @title Check linkage colors
##' @param geneVec A vector containing condidate genes used for Circos plot. The length of one is also allowed.
##' @param linkColVec A vector of colors used for linkages.
##' @param allAnno The annotation matrix. The 1st is the KEGG gene IDs, which has the same format with "inputft". The 2nd is the symbol name shown in the Circos plot. The 3rd column is the chromosome name. The 4th and 5th columns are the start and end positions. The "allAnno" may contains genes that are absent in the "inputft".
##' @return A list containing warning message.
##' @examples
##' data(geneAnno)
##' ## no annotation and containing NA
##' g1 <- c('hsa:111111', 'hsa:498', 'hsa:506', 'hsa:509', 'hsa:516', 'hsa:517')
##' g1Col <- c('blue', 'dblue_a3', 'dyellow_a3', 'dred_a3', NA, NA)
##' g1 <- c(NA, NA)
##' g1 <- c('hsa:0000', 'hsa:111111')
##' g1 <- c('hsa:111111', NA, NA)
##' g1 <- c('hsa:111111', 'hsa:498', NA, NA)
##' checkColList <- CheckLinkCol(g1, g1Col, geneAnno)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
CheckLinkCol <- function(geneVec, linkColVec, allAnno) {

  ## initiate
  checkGeneVec <- character(0)
  checkLinkCol <- character(0)
  checkGeneSym <- character(0)
  wm <- character(0)
  ## warning message template
  wmTemp <- 'No legal colours or annotated genes are found.\n'

  hasLinkCol <- linkColVec[!is.na(linkColVec)]
  geneVec <- geneVec[!is.na(linkColVec)]

  if (length(hasLinkCol) == 0) {
    ## no legal linkage colours.
    wm <- c(wm, wmTemp)
  } else {
    hasLogic <- geneVec %in% allAnno[, 1]
    hasAnnoGeneVec <- geneVec[hasLogic]
    if (length(hasAnnoGeneVec) == 0) {
      ## no annotated genes
      wm <- c(wm, wmTemp)
    } else {
      checkGeneVec <- hasAnnoGeneVec
      checkLinkCol <- hasLinkCol[hasLogic]
      checkGeneSym <- allAnno[match(checkGeneVec, allAnno[, 1]), 2]
    }
  }

  reList <- list(checkGeneVec = checkGeneVec,
                 checkLinkCol = checkLinkCol,
                 checkGeneSym = checkGeneSym,
                 wm = wm)

  return(reList)
  
}


##' Selection and Annotation of from-to matrix
##'
##' Annotation the "from-to matrix"(network data). The main purpose is: 1. select ft matrix containing certain gene list; 2. transfer geneID (hsa:1) to gene symbol (A1BG).
##' @title Selection and Annotation of ft matrix
##' @inheritParams writeCircos
##' @return An annotated ft matrix
##' @examples
##' data(atpft)
##' atpft <- atpft[, c(1, 3, 5:6)]
##' data(geneAnno)
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' newft <- Annoft(f1genes, atpft, geneAnno)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
Annoft <- function(geneVec, inputft, allAnno) {

  ## select ft matrix containing the 'geneVec'
  ftLogic <- (inputft[, 1] %in% geneVec) | (inputft[, 2] %in% geneVec)
  inputft <- inputft[ftLogic, ]

  ## anno from and to genes
  annoFrom <- allAnno[match(inputft[, 1], allAnno[, 1]), 2]
  annoTo <- allAnno[match(inputft[, 2], allAnno[, 1]), 2]
  inputft <- cbind(annoFrom, annoTo, inputft[, -1:-2])

  ## remove NA
  annoLogic <- (is.na(inputft[, 1])) | (is.na(inputft[, 2]))
  inputft <- inputft[!annoLogic, ]

  return(inputft)
}

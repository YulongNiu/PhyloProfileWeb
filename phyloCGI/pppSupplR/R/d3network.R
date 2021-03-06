##' Plot network using the javascript D3 library
##'
##' D3 library is used for network plot.
##' @title plot D3 network
##' @param d3Annoft A ft matrix used for d3 plot. 1st column: from nodes. 2nd column: to nodes. 3rd column: linkage strength, for example the Jaccard similarity ranging from 0 to 1. 4th column: linkage colours.
##'
##' A normalized linkage strength is used. The nodes are grouped by different color, and the latter gene will cover the former one. The size of a node is determined by the number of its partners.
##' @return d3Transft returns a list, the first element is links dataframe and the second element is the nodes dataframe.
##' PlotNetworkD3 returns an object of D3 network.
##' writeD3network returns nothing.
##' @examples
##' filePath <- system.file("extdata", "linkage.csv", package = "pppSupplR")
##' atpft <- read.csv(filePath, stringsAsFactor = FALSE)
##' atpft <- atpft[, c(4, 8, 9, 11)]
##' newd3ft <- d3Transft(atpft)
##' d3Obj <- d3PlotNet(newd3ft)
##' \dontrun{
##' writed3Net(d3Obj, 'testd3.html')
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname plotd3
##' @export
##'
d3Transft <- function(d3Annoft) {
  ## entire gene set
  fSet <- unique(d3Annoft[, 1])
  tSet <- unique(d3Annoft[, 2])
  entireSet <- unique(c(fSet, tSet))
  entireLen <- length(entireSet)

  ## construct link data frame
  ## from and to index
  fromIdx <- match(d3Annoft[, 1], entireSet)
  toIdx <- match(d3Annoft[, 2], entireSet)
  linkStrength <- as.numeric(d3Annoft[, 3])

  d3Links <- data.frame(fromIdx = fromIdx - 1,
                        toIdx = toIdx - 1,
                        linkStrength = linkStrength)

  ## construct node data frame
  colorVec <- character(entireLen)
  for (i in seq_along(fSet)) {
    eachMat <- d3Annoft[d3Annoft[, 1] == fSet[i], , drop = FALSE]
    eachNode <- c(eachMat[, 2], fSet[i])
    eachLogic <- entireSet %in% eachNode
    colorVec[eachLogic] <- ifelse(colorVec[eachLogic] == '',
                                  eachMat[1, 4],
                                  colorVec[eachLogic])
  }

  ## node size is denoted by the number of its partners
  sizeVec <- numeric(entireLen)
  for (i in 1:entireLen) {
    interLogic <- (d3Annoft[, 1] == entireSet[i]) | (d3Annoft[, 2] == entireSet[i])
    sizeVec[i] <- sum(interLogic)
  }

  d3Nodes <- data.frame(entireSet = entireSet,
                        sizeVec = sizeVec,
                        colorVec = colorVec)
  ## d3 ft list
  d3ftList <- list(d3Links = d3Links, d3Nodes = d3Nodes)

  return(d3ftList)
}


##' @param d3ft The d3 ft list from d3Transft().
##' @param linkStrengthBase The base value of linkage strength.
##' @param ... Additional paramters from forceNetwork().
##' @rdname plotd3
##' @importFrom networkD3 forceNetwork JS
##' @export
##'
d3PlotNet <- function(d3ft,
                      linkStrengthBase = 5,
                      ...) {

  d3Links <- d3ft[[1]]
  d3Nodes <- d3ft[[2]]

  ## linkage strength
  d3Links[, 3] <- linkStrengthBase * d3Links[, 3]
  ## d3Links[, 3] <- linkStrengthBase * Normalize(d3Links[, 3], addSmall = 1e-5)
  ## d3Nodes[, 3] <- nodeSizeBase * Normalize(d3Nodes[, 3], addSmall = 1e-5)

  d3netObj <- forceNetwork(Links = d3Links, Nodes = d3Nodes,
                           Source = 'fromIdx', Target = 'toIdx',
                           Value = 'linkStrength', NodeID = 'entireSet',
                           Group = 'colorVec', Nodesize = 'sizeVec',
                           opacity = 1, linkColour = '#dfdfdf',
                           charge = -30, bounded = TRUE,
                           height = 1000, width = 1000,
                           zoom = TRUE, fontSize = 18,
                           colourScale = JS("e = function colourScale(colour) { return colour;}"), ...)

  return(d3netObj)
}

##' @param d3NetObj A d3 network object from PlotNetworkD3().
##' @param fileName The file name to store the html.
##' @param savePath The save path of html.
##' @importFrom networkD3 saveNetwork
##' @rdname plotd3
##' @export
##'
writed3Net <- function(d3NetObj,
                       fileName = 'networkd3.html',
                       savePath = './') {

  saveFullPath <- paste0(savePath, fileName)

  saveNetwork(d3NetObj, file = saveFullPath, selfcontained = FALSE)
}


##' Extract d3 network elements
##'
##' Retrieve the d3 network elements from a html file.
##' @title Extract d3 network
##' @param d3htmlPath An absolute path of html, for example '/home/User1/Download/index.html'.
##' @return A string vector containing the d3 network html elements.
##' @examples
##' filePath <- system.file("extdata", "d3net.html", package = "pppSupplR")
##' d3Ele <- d3ExtractNetEle(filePath)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @export
##'
d3ExtractNetEle <- function(d3htmlPath) {

  ## file path and read in
  filePath <- paste0('file://', d3htmlPath)
  entireD3HTML <- getURL(filePath)

  ## d3 content
  ## first <div> to fourth </script>
  divReg <- gregexpr('<div.*?</div>', entireD3HTML)
  scriptReg <- gregexpr('<script.*?</script>', entireD3HTML)
  d3Val <- substring(entireD3HTML,
                     divReg[[1]],
                     scriptReg[[1]][5] - 2)

  return(d3Val)
}


##' Normalization of the numeric vector
##'
##' The normalization formular is (x - min(X))/(max(X) - min(X)). "NA"s are removed. If the max(x) is equal to min(x), asign all the elements to 0.5.
##' @title normalization
##' @param inputVec The input raw vector.
##' @param addSmall A small value used to add on the normalized data.
##' @return normalized vector
##' @examples
##' testVec <- sample(1:100, 10)
##' Normalize(testVec)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
Normalize <- function(inputVec, addSmall = 0) {

  maxVal <- max(inputVec, na.rm = TRUE)
  minVal <- min(inputVec, na.rm = TRUE)

  if (maxVal == minVal) {
    normVec <- rep(0.5, length(inputVec))
  } else {
    normVec <- (inputVec - minVal)/(maxVal - minVal)
  }

  ## add small value
  normVec <- normVec + addSmall

  return(normVec)
}


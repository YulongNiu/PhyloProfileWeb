## May be deleted in future.

##' @param linkColVec A vector of colors used for linkages, names, and heatmaps.
##' @param geneSymbolVec A vector containing condidate genes symbol.
##' @rdname circosConf
##' @return A formated rules vector.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
cRules <- function(geneSymbolVec, linkColVec) {
  ## template 
  ruleTemplate <- '<rule>\r\nimportance = %importValue%\r\ncondition = var(value) eq "%geneName%"\r\ncolor = %linkCol%\r\n</rule>\r\n'

  geneLen <- length(geneSymbolVec)
  importValVec <- seq(90, 95 - geneLen * 5, by = -5)
  ruleList <- list()

  for(i in 1:geneLen) {
    geneRule <- sub('%importValue%', importValVec[i], ruleTemplate)
    geneRule <- sub('%geneName%', geneSymbolVec[i], geneRule)
    geneRule <- sub('%linkCol%', linkColVec[i], geneRule)
    ruleList[[i]] <- geneRule
  }

  ruleVec <- paste(unlist(ruleList), collapse = '')

  rulesVec <- paste0('\r\n<rules>\r\n', ruleVec, '</rules>\r\n')
  return(rulesVec)
}



##' Generate configures for Circos plot
##'
##' "labelConf" generates gene labels for the Circos configures ("names and locations of genes"). "cHeatmap" generates heat map for the Circos configures ("heatmap"(correlation value) for each linkage). "cFreq" generates frequency of presence in three domains. "cLinks" generates linkages for each input "geneVec". "cPhyloLink" generates the whole file automatically.
##' Please NOTE this temporal configure code is a piece of shit. Please check json format and re-organize this shit.
##' @title Generate Circos configure
##' @param savePath The folder of label files.
##' @inheritParams cRules
##' @return A formated label vector.
##' @rdname circosConf
##' @examples
##' f1 <- c('ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5O', 'ATP5D', 'ATP5E')
##' f1Col <- c('dblue_a3', 'dyellow_a3', 'dred_a3', 'dgreen_a3', 'dpurple_a3', 'dorange_a3')
##' labelConf <- cLabelGenes('phylo/', f1, f1Col)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
cLabelGenes <- function(savePath, geneSymbolVec, linkColVec) {

  ## rules
  ruleConf <- cRules(geneSymbolVec, linkColVec)

  ## temple
  labelTemp <- '\r\n<plot>\r\ntype = text\r\ncolor = black\r\nfile = %labelFile%labelGene.txt\r\nr0 = 1.04r\r\nr1 = 1.04r+500p\r\nshow_links = yes\r\nlink_dims = 2p,2p,2p,2p,2p\r\nlink_thickness = 1p\r\nlink_color = red\r\npadding = 2p\r\nrpadding = 2p\r\nlabel_snuggle = yes\r\nmax_snuggle_distance = 1r\r\nsnuggle_tolerance = 0.25r\r\nsnuggle_sampling = 2\r\nsnuggle_link_overlap_test = yes\r\nsnuggle_link_overlap_tolerance = 2p\r\nsnuggle_refine = yes\r\nlabel_rotate = yes\r\nlabel_size = 13p\r\nlabel_font = condensed\r\n%rules%\r\n</plot>\r\n'

  labelVec <- sub('%labelFile%', savePath, labelTemp)
  labelVec <- sub('%rules%', ruleConf, labelVec)
  
  return(labelVec)
}



##' @param geneVec A vector containing condidate genes which should be the same as the input of function "writeCircos". CAUTION: The element in the "geneVec" will be used in the files' name!!! An example is c('hsa:516', 'hsa:517')
##' @inheritParams cLabelGenes
##' @seealso writeCircos
##' @return A formated heatmap vector.
##' @rdname circosConf
##' @examples
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' f1Col <- c('dblue_a3', 'dyellow_a3', 'dred_a3', 'dgreen_a3', 'dpurple_a3', 'dorange_a3')
##' heatConf <- cHeatmap('phylo/', f1genes, f1Col)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
cHeatmap <- function(savePath, geneVec, linkColVec) {

  geneLen <- length(geneVec)
  
  ## heattemp
  heatTemp <- '\r\n<plot>\r\ntype = scatter\r\nstroke_thickness = 1\r\nfile = %heatFile%\r\ncolor = %linkCol%\r\nstroke_color = black\r\nglyph = circle\r\nglyph_size = 15\r\nmax = 1\r\nmin = 0\r\nr1 = 0.96r\r\nr0 = 0.85r\r\n'

  ## background temp
  bgTemp <- '\r\n<backgrounds>\r\n<background>\r\ncolor = vvlgrey\r\ny0 = 0\r\ny1 = 1.01\r\n</background>\r\n</backgrounds>\r\n<axes>\r\n<axis>\r\ncolor = lgrey\r\nthickness = 1\r\nspacing = 0.2r\r\ny0 = 0.0\r\ny1 = 1.01\r\n</axis>\r\n</axes>\r\n</plot>\r\n'
  
  fileNames <- paste0(savePath, geneVec, 'heatmap.txt')
  heatList <- list()

  for(i in 1:geneLen) {
    if (i == 1) {
      heatEach <- paste0(heatTemp, bgTemp)
    } else {
      heatEach <- paste0(heatTemp, '</plot>\r\n')
    }
    
    heatEach <- sub('%heatFile%', fileNames[i], heatEach)
    heatEach <- sub('%linkCol%', linkColVec[i], heatEach)

    heatList[[i]] <- heatEach
  }

  heatVec <- paste(unlist(heatList), collapse = '')
  return(heatVec)
  
}


##' @inheritParams cHeatmap
##' @return A formated frequency
##' @rdname circosConf
##' @examples
##' freqConf <- cFreq('phylo/')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
cFreq <- function(savePath) {

  freqTemp <- '\r\n<plot>\r\ntype = scatter\r\nstroke_thickness = 1\r\nfile = %filePath%freqArc.txt\r\ncolor = vlgrey\r\nstroke_color = black\r\nmin = 0\r\nmax = 1\r\nr0 = 0.80r\r\nr1 = 0.83r\r\n<backgrounds>\r\n<background>\r\ncolor = vvlgrey\r\ny0 = 0\r\ny1 = 1\r\n</background>\r\n</backgrounds>\r\n</plot>\r\n<plot>\r\ntype = scatter\r\nstroke_thickness = 1\r\nfile = %filePath%freqBac.txt\r\ncolor = lgrey\r\nstroke_color = black\r\nmin = 0\r\nmax = 1\r\nr1 = 0.79r\r\nr0 = 0.76r\r\n<backgrounds>\r\n<background>\r\ncolor = vlgrey\r\ny0 = 0\r\ny1 = 1\r\n</background>\r\n</backgrounds>\r\n</plot>\r\n<plot>\r\ntype = scatter\r\nstroke_thickness = 1\r\nfile = %filePath%freqEu.txt\r\ncolor = grey\r\nstroke_color = black\r\nmin  = 0\r\nmax = 1\r\nr1 = 0.75r\r\nr0 = 0.72r\r\n<backgrounds>\r\n<background>\r\ncolor = lgrey\r\ny0 = 0\r\ny1 = 1\r\n</background>\r\n</backgrounds>\r\n</plot>'

  freqVec <- gsub('%filePath%', savePath, freqTemp)
  return(freqVec)
}


##' @inheritParams cHeatmap
##' @return A formated links
##' @rdname circosConf
##' @examples
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' f1Col <- c('dblue_a3', 'dyellow_a3', 'dred_a3', 'dgreen_a3', 'dpurple_a3', 'dorange_a3')
##' linkConf <- cLinks('phylo/', f1genes, f1Col)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
cLinks <- function(savePath, geneVec, linkColVec) {

  geneLen <- length(geneVec)
  fileNames <- paste0(savePath, geneVec, '.txt')

  linkTemp <- '\r\n<link>\r\nz = 50\r\ncolor = %linkCol%\r\nfile = %linkFile%\r\nbezier_radius_purity = 0.2\r\ncrest = 1\r\n</link>'
  linkList <- list()

  for(i in 1:geneLen) {
    eachLink <- sub('%linkCol%', linkColVec[i], linkTemp)
    eachLink <- sub('%linkFile%', fileNames[i], eachLink)
    linkList[[i]] <- eachLink
  }
  linkVec <- paste(unlist(linkList), collapse = '')

  linksTemp <- '\r\n<links>\r\nz = 0\r\nradius = 0.7r\r\ncrest = 0.5\r\nbezier_radius = 0.5r\r\nbezier_radius_purity = 0.75\r\n%linkData%\r\n</links>\r\n'
  linksVec <- sub('%linkData%', linkVec, linksTemp)

  return(linksVec)
}




##' @inheritParams cHeatmap
##' @inheritParams cLabelGenes
##' @param organism KEGG species ID, for example "hsa" for human.
##' @return A formated Circos configure vecor of predicted linkages.
##' @rdname circosConf
##' @examples
##' f1 <- c('ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5O', 'ATP5D', 'ATP5E')
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' f1Col <- c('dblue_a4', 'dyellow_a4', 'dred_a4', 'dgreen_a4', 'dpurple_a4', 'dorange_a4')
##' wholeConf <- cPhyloLink(savePath = 'phylo/', geneVec = f1,
##' geneSymbolVec = f1, linkColVec = f1Col, organism = 'hsa')
##' wholeConf <- cPhyloLink(savePath = 'phylo/', geneVec = f1genes,
##' geneSymbolVec = f1, linkColVec = f1Col, organism = 'hsa')
##' 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
cPhyloLink <- function(savePath, geneVec, geneSymbolVec, linkColVec, organism) {

  phyloTemp <- '<<include etc/colors_fonts_patterns.conf>>\r\n<<include ideogram.conf>>\r\n<<include ticks.conf>>\r\n<image>\r\n<<include etc/image.conf>>\r\nradius* = 3000p\r\n</image>\r\nchromosomes_units = 1000000\r\nchromosomes_display_default = yes\r\nkaryotype = karyotype.%organismData%.txt\r\n<plots>%labelData%\r\n%heatmapData%\r\n%freqData%\r\n</plots>%linkData%\r\n<<include etc/housekeeping.conf>>\r\n'

  labelConf <- cLabelGenes(savePath, geneSymbolVec, linkColVec)
  heatmapConf <- cHeatmap(savePath, geneVec, linkColVec)
  freqConf <- cFreq(savePath)
  linkConf <- cLinks(savePath, geneVec, linkColVec)

  phyloVec <- sub('%labelData%', labelConf, phyloTemp)
  phyloVec <- sub('%heatmapData%', heatmapConf, phyloVec)
  phyloVec <- sub('%freqData%', freqConf, phyloVec)
  phyloVec <- sub('%linkData%', linkConf, phyloVec)
  phyloVec <- sub('%organismData%', organism, phyloVec)

  return(phyloVec)
}



##' @inheritParams cPhyloLink
##' @param storePath The store path of Circos configure file.
##' @return NULL
##' @rdname circosConf
##' @examples
##' f1 <- c('ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5O', 'ATP5D', 'ATP5E')
##' f1genes <- c('hsa:498', 'hsa:506', 'hsa:509', 'hsa:539', 'hsa:513', 'hsa:514')
##' f1Col <- c('dblue_a4', 'dyellow_a4', 'dred_a4', 'dgreen_a4', 'dpurple_a4', 'dorange_a4')
##' \dontrun{
##' writeConf(savePath = 'phylo/', geneVec = f1genes,
##' geneSymbolVec = f1, linkColVec = f1Col,
##' organism = 'hsa', storePath = 'tmp1/')
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
writeConf <- function(savePath, geneVec, geneSymbolVec, linkColVec, organism, storePath) {

  confVec <- cPhyloLink(savePath, geneVec, geneSymbolVec, linkColVec, organism)
  fileConn <- file(paste0(storePath, 'circosConf.conf'))
  writeLines(confVec, fileConn)
  close(fileConn)

  return(NULL)
}



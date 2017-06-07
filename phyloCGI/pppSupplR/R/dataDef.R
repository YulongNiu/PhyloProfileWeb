##' Whole genome predicted linkages of human F1Fo ATP synthase subunits.
##'
##' A character matrix:
##' the 1st column and 2ed colum are gene ID;
##' the 3rd clumn is the Jaccard similarity ranging from 0 to 1.
##'
##' @docType data
##' @name atpft
##' @format A character matrix
##' @references Unpublished data from Yulong Niu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##'
NULL


##' Annotatio of human genes.
##'
##' A character matrix:
##' the 1st column is the gene ID which should be the same format in "atpft";
##' the 2nd column is the chromosome name, like "1", "2", and "MT";
##' the 3rd and 4th columns are the start and end position, respectively.
##'
##' @docType data
##' @name geneAnno
##' @format A character matrix
##' @references Unpublished data from Yulong Niu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##'
NULL


##' Phylogeny of 972 KEGG species
##'
##' A character matrix:
##' the 1st column is KEGG species IDs;
##' the 2nd column is the corresponding phylogeny;
##'
##' @docType data
##' @name phyloSpe
##' @format A character matrix
##' @references Unpublished data from Yulong Niu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##'
NULL


##' Phylogenic profiles
##'
##' A character matrix:
##' the column names are the KEGG species IDs;
##' the row names are the KEGG gene IDs;
##'
##' @docType data
##' @name wholeProfile
##' @format A character matrix
##' @references Unpublished data from Yulong Niu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##'
NULL

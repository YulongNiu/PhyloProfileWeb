# User's guide #

<!-- content start -->

**Table of Contents**

- [1. Input](#1-input)
    - [1.1 Linkage number](#11-linkage-number)
    - [1.2 Plot parameters](#12-plot-parameters)
    - [1.3 Candidate gene list](#13-candidate-gene-list) 
- [2. Output](#2-output)

<!-- content end -->


## 1. Input ##

### 1.1 Linkage number ###

The threshold of linkages number for each gene could be randomly chose between `1` and `500`(even a float number like `5.4`). A smaller vaule means a more stringent threshold of predicted linkages.

If you want to get less but more reliable linkages, please try a small value, for example `20`. Please be aware, a small value may be at a risk of missing some potential and novel linkages which do not have very strong evidence.

### 1.2 Plot parameters ###

Please try default plot parameters at first. If you are not satisfied with output figures, please adjust the parameters and run it again.

### 1.3 Candidate gene list ###

Only the [`txt`](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.txt) (with the separator of `"Tab"`) or the [`csv`](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.csv) (with the separator of `","`) format is allowed for the input files. Please do not try to upload other types of files (like `test.html`). Otherwise, an error web-page may return.

The number of input genes is supposed to be at least two. Otherwise, an error message will return. There are several reasons. The most important one is our algorithm is based on the whole genome, and the prediction results are more reliable when taking into account a cluster of genes like genes from a complex or a biological pathway. Our team is now working on an advanced algorithms based on the statistical distribution of linkages for each gene. Several clues indicate the new algorithm have a better performace in reducing the false positive rate, which haunts all the existing methods for protein linkage prediction.

* 1st column

1st column contains the input genes. The "GeneID" is used as the unique gene name for each gene, for example `hsa:513` denotes the human gene "ATP5D" (ATP synthase, H+ transporting, mitochondrial F1 complex, delta subunit (EC:3.6.1.14)). A total of 20127 [human whole genomic genes](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/wholeGenomeGenes.csv), **please refer this to set your own input genes**. Other formats of gene names may be supported in the future.

* 2nd column

2nd column contains the colours for each input genes. Colour modes like the hex triplet and colour names are both supported. For example, the "red" and `#ff0000` are legal input colours. If you do not want to show the colour for a gene, an easy way is to set the colour as `white` or `#FFFFFF`.

* 3rd column

3rd column defines linkages colours. The colour format is different from the ones in the 2nd column. The legal input colours are described as the [Circos configuration](http://circos.ca/documentation/tutorials/configuration/colors/). Seven colours, red, orange, yellow, green, blue, purple, and grey are predefined with prefix to show the light (`l`, `vl`, `vvl`) and dard (`d`, `vd`, `vvd`) verions. Five degrees of [transparent colous](http://circos.ca/tutorials/lessons/recipes/transparent_links/) are allowed, `_a1` have a 17% transparency, `_a2` 33%, `_a3` 50%, `_a4` 67% and `_a5` 83%. For example, `vlred_a4` is the very light red colour with 67% transparency. You may want to hidden the linkaged for certain genes, just use `NA` (not "na", not "white", not "Na") as the link color. But, at least one genes should have proper linkage colour.

An internal [gene annotation file](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/geneAnno.csv) is used for Circos plot. **Any genes not in this list will be neglected**. The time of plotting the Circos figure increases dramatically along with increment of predicted linkages. As the web-server connection will be closed in 15 minutes, the input candidate gene number for Circos plot  we recommend is no more than `7`.

* 4th column

4th column allows pre-defined gene names (optional). As the GeneID (*e.g.*, `hsa:513`) is hard to get its meaning, the 4th column could provide a way to transfer the GeneID to other formats ([Example](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.txt)). The pre-defined gene names will be used inboth phylogenetic profile and correlation figures. This column is optional; if you want to keep the GeneID in figures, please leave it as blank (no column names and no values).

## 2. Output ##

We provide two ways to illustrate the prediction results. A webpage merged the all the figures and tables will be returned at first. At the bottom of the webpage, a linkage will guide you to download all the results into your local computer.

* Webpage

The webpage contains visualization of phylogenetic profiles, correlation matrix, and list of prediction linkages. A grey block will be used to represent the negative correlation value in the correlation matrix plot.

* Download folder

The downloaded folder includes published quality of figures mentioned above in both "pdf" and "jpg" formats. The correlation matrix and prediction linkages is stored in a "csv" file that can be easily used for other analysis and validation experiments. Moreover, the same webage is stored in a "html" file used for a overview of the results.

Please download and save your results (at the end of the return webpage), the result files will be kept on the web-server only `2` hours.








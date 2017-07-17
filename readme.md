# User's guide #

<!-- content start -->

**Table of Contents**

- [1. Input](#1-input)
    - [1.1 Algorithm parameters](#11-algorithm-parameters)
    - [1.2 Plot parameters](#12-plot-parameters)
    - [1.3 Candidate gene list](#13-candidate-gene-list) 
    - [1.4 Annotation files](#14-annotation-files) 
- [2. Output](#2-output)

<!-- content end -->


## 1. Input ##

### 1.1 Algorithm parameters ###

* Linkages number

The threshold of linkages number for each gene could be randomly chose between `1` and `500`(even a float number like `5.4`). A smaller vaule means a more stringent threshold of predicted linkages.

If you want to get less but more reliable linkages, please try a small value, for example `20`. Please be aware, a small value may be at a risk of missing some potential and novel linkages which do not have very strong evidence.

* BLAST Evalue

Three BLAST E-value threshold, 0.001, 0.0005, 0.0001, were provided. This threshold was used in choosing homologous proteins. For a pair proteins (the referenced protein and test protein), if the blast E-value that is smaller than the threshold, the homologs is through to be present; otherwise, the homologs is absent. Correspondingly, in the phylogenetic profile, "1" and "0" are used to denote the presence and absence.

* Organism

*Homo sapiens* (human) and *Arabidopsis thaliana* (thale cress) are now supported. 

### 1.2 Plot parameters ###

Please try default plot parameters at first. If you are not satisfied with output figures, please adjust the parameters and run it again.

### 1.3 Candidate gene list ###

Only the [`txt`](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.txt) (with the separator of `"Tab"`) or the [`csv`](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.csv) (with the separator of `","`) format is allowed for the input files. Please do not try to upload other types of files (like `test.html`). Otherwise, an error webpage may return.

The number of input genes is supposed to be at least two. Otherwise, an error message will return. There are several reasons. The important one is our algorithm is based on the whole genome, and the prediction results are more reliable when taking into account a cluster of genes like genes from a complex or a biological pathway.

* 1st column

1st column contains the input genes. The "GeneID" is used as the unique gene name for each gene, for example `hsa:513` denotes the human gene "ATP5D" (ATP synthase, H+ transporting, mitochondrial F1 complex, delta subunit (EC:3.6.1.14)). **Please refer the "[Input gene list](#14-annotation-files)" to set your own input genes**. Other formats of gene names may be supported in the future.

* 2nd column

2nd column contains the colours for each input genes. Colour modes like the hex triplet and colour names are both supported. For example, the "red" and `#ff0000` are legal input colours. If you do not want to show the colour for a gene, an easy way is to set the colour as `white` or `#FFFFFF`.

* 3rd column

3rd column defines linkages colours, which uses the same colour format as the 2nd column. You may want to hidden the linkaged for certain genes, just use `NA` (not "na", not "white", not "Na") as the link color. But, at least one genes should have proper linkage colour.**An internal [Circos gene list](#14-annotation-files) is used for Circos plot, and any genes not in this list will be neglected**. The D3 interactive network is also constructed based on the 3rd column data. It is designed as an alternative way to show the linkages network.

* 4th column

4th column allows pre-defined gene names (optional). As the GeneID (*e.g.*, `hsa:513`) is hard to get its meaning, the 4th column could provide a way to transfer the GeneID to other formats ([Example](http://bioinfor.scu.edu.cn/phyloprofile/Exampledata/phylopred_fatp1Link/atpSubOne.txt)). The pre-defined gene names will be used inboth phylogenetic profile and correlation figures. This column is optional; if you want to keep the GeneID in figures, please leave it as blank (no column names and no values).

### 1.4 Annotation files ###

|Organism|Input gene list (1st column)|Circos gene list|
|:------:|:--------------------------:|:----------:|
|*Homo sapiens*|20127 [Download](http://bioinfor.scu.edu.cn/phyloprofile/AnnoData/hsa_wholeGenomeAnno.csv)|19077 [Download](http://bioinfor.scu.edu.cn/phyloprofile/AnnoData/hsa_geneAnno.csv)|
|*Arabidopsis thaliana*|27396 [Download](http://bioinfor.scu.edu.cn/phyloprofile/AnnoData/ath_wholeGenomeAnno.csv)|27394 [Download](http://bioinfor.scu.edu.cn/phyloprofile/AnnoData/ath_geneAnno.csv)|

## 2. Output ##

We provide two ways to illustrate the prediction results. A webpage merged the all the figures and tables will be returned at first. At the bottom of the webpage, a linkage will guide you to download all the results into your local computer.

* Webpage

The webpage contains visualization of phylogenetic profiles, correlation matrix, Circos plot, D3 interactive network, and list of prediction linkages. A grey block will be used to represent the negative correlation value in the correlation matrix plot.

* Download folder

The downloaded folder includes published quality of figures mentioned above in both "pdf" and "jpg" formats. The correlation matrix and prediction linkages are stored in a "csv" file that can be easily used for other analysis and validation experiments. Moreover, the same webage is stored in a "html" file used for a overview of the results.

Please download and save your results (at the end of the return webpage), the result files will be kept on the web-server only `2` hours.








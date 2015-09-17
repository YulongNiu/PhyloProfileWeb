# User's guide #

<!-- content start -->

**Table of Contents**

- [1. Input](#1-input)
    - [1.1 Linkage number](#11-linkage-number)
    - [1.2 Plot parameters](#12-plot-parameters)
- [2. Output](#2-output) 
<!-- content end -->


## 1. Input ##

### 1.1 Linkage number ###

The threshold of linkages number for each gene could be randomly chose between `1` and `500`(even a float number like `5.4`). A smaller vaule means a more stringent threshold of predicted linkages.

If you want to get less but more reliable linkages, please try a small value, for example `20`. Please be aware, a small value may be at a risk of missing some potential and novel linkages which do not have very strong evidence.

### 1.2 Plot parameters ###

Please try default plot parameters at first. If you are not satisfied with output figures, please adjust the parameters and run it again.

## 2. Output ##

We provide two ways to illustrate the prediction results. A webpage merged the all the figures and tables will be returned at first. At the bottom of the webpage, a linkage will guide you to download all the results into your local computer.

* Webpage

The webpage contains visualization of phylogenetic profiles, correlation matrix, and list of prediction linkages. A grey block will be used to represent the negative correlation value in the correlation matrix plot.

* Download folder

The downloaded folder includes published quality of figures mentioned above in both "pdf" and "jpg" formats. The correlation matrix and prediction linkages is stored in a "csv" file that can be easily used for other analysis and validation experiments. Moreover, the same webage is stored in a "html" file used for a overview of the results.

Please download and save your results (at the end of the return webpage), the result files will be kept on the web-server only `2` hours.








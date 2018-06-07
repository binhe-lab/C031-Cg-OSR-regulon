---
title: Suggestions for starting the analysis
author: Bin He
date: 2018-06-06
---

# Useful resources
- [Data from Roetzer et al 2011](https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2915/)
    There are three types of data you can download: the raw data, the processed data and the R expression set. The last one is the most "processed" and allows you to work directly with the formatted data. You can start with this to make the learning curve less steep. Feel free to explore the raw and processed dataset.
- [working with R expression set](https://bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf)
    This shows you how to work with the eSet data in R
- [Tutorial on microarray analysis using R](http://legacydirs.umiacs.umd.edu/~hcorrada/CFG/lectures/lect21_bioc/bioconductor.pdf)
- [Limma user guide](http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
    Limma is a super useful package for analyzing microarray and RNA-seq data. The user guide is not only a manual, but also contains a wealth of information about the theory. There are references in it that explains the details of the algorithms.
- How to install Bioconductor packages
    ```r
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma")
    ```

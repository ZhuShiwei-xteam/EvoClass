# EvoClass

#### Contents

* [Overview](#overview)
* [Repo Contents](#repo-contents)
* [Instructions for Use](#instructions-for-use)
  * [Obtain data](#obtain-data)
  * [Reproduce analysis](#reproduce-analysis)
* [Test Environment](#test-environment)

## Overview

Purpose of this repository is to share analysis procedure, data and help readers or reviewers to know more detail of this work, reproduce or make use of results they are interested in.


## Repo Contents

* [code](./code): tidy R functions and R script.
* [data](./data): original and preprocessed data used for analysis and share.
* [result](./result): important middle results and final results of manuscript, most of them are in form of `.RData`, which can be easily loaded and operated by R. 

## Instructions for Use

### Obtain data

For readers who want to obtain raw/result data, locate data file, then download it with one of following ways:

* In Github, download file by clicking either `Download` button or `Raw` button at corresponding data page

* Use linux command `wget` or `curl`, fo example, you can download PCNA gene sets by

  `wget https://github.com/ZhuShiwei-xteam/EvoClass/blob/master/data/PCNA_Signature.txt`

Or you can download whole respository with one of following ways:

* Clone this repository with `git clone https://github.com/ZhuShiwei-xteam/EvoClass.git`
* Download whole respository by clicking `Download` button at top right of url page <https://github.com/ZhuShiwei-xteam/EvoClass>

### Reproduce analysis

For readers who want to reproduce analysis shown in manuscript, please [install R](https://cran.r-project.org) in your computer.

## Test Environment

* System: __Linux__

* Software: __R v3.6.2__

* R packages:
  * ggsci
  * [tidyverse](https://www.tidyverse.org/)
  * oncosign
  * ggpubr
  * pheatmap
  * mclust
  * survival
  * survminer
  * powerSurvEpi
  * ggplot2
  * AnnotationDbi
  * org.Hs.eg.db
  * DESeq2
  * dplyr
  * data.table
  * clusterProfiler
  * sequenza

These R packages are easily searched by internet and have no strict version requirements to reproduce the analyses.


* Other Tools: __PyClone v0.13.0; GSEA Desktop v4.1.0 requires Java 8__


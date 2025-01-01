# **Gene Expression Analysis and Visualisation of a RNA-seq dataset**

## 1. **Project Objective**

In this project, a processed RNA-seq dataset with quantitated gene expression data from an RNA-seq experiment was analyzed to identify the significantly expressed genes. 

## 2. **GitHub repository folders**

* Code \- contains the R script of the code for data preprocessing and DEA.  
* Data \- contains the RNA-seq dataset.  
* Results \- contains the data and results obtained after pre-processing, differential expression analysis (DEA) and visualization.


## 3. **Requirements**

   The following R libraries were used:
* dplyr    
* ggplot2    
  
These can be installed by running:  
install.packages(c(“dplyr”, “ggplot2”))

## 4. **Methodology**
* Download the following RNA-seq dataset: https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt   
* Read the data using **read.table()** function.
* Generate a volcano plot for the differentially expressed genes (DEGs) using the **ggplot()** function.   
* Determine the up-regulated genes (genes with log2FC >1 and p-value <0.01).
* Determine the down-regulated genes (genes with log2FC <-1 and p-value <0.01) 
* The functions of the top 5 significantly up-regulated and down-regulated genes were obtained using the [GeneCards](https://www.genecards.org/) database.

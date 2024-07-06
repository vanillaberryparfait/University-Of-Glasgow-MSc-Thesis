# MSc Thesis Project as Carried Out At The University Of Glasgow 2022-2023

This repository contains the scripts and data used in my MSc thesis project. In this project , I investigated the efficacy of combining the BTK inhibitor ibrutinib and the mTOR inhibitor AZD8055 in treating chronic lymphocytic leukemia (CLL). By analyzing differentially expressed genes in CLL patients under various conditions, I explored synergistic effects, enriched pathways, and potential links between CLL and neurodegenerative disorders

## Table of Contents

- [Data](#data)
- [Installation](#installation)
- [Workflow](#workflow)
- [Overview](#overview)

## Data

The data used in this study was a gene count data obtained via Next Generation
Sequencing (NGS) which was done at Novogene. This dataset was provided to us by Dr.
Alison Michie from the Paul O'Gorman Leukaemia Research Centre. The dataset
encompassed gene expression information from five male patients from the UK, each
uniquely identified by patient identifiers, namely CLL_125, CLL_173, CLL_175,
CLL_186, and CLL_203. All of these patients were stated to be older adults i.e., above
50 years of age. The patients were classified into different Binet stages based on severity.
Patient CLL_203, CLL_173 and CLL_125 were observed to have mutations and
cytogenetic abnormalities. Patient CLL_203 did not show 11q, 17p or p53 mutations,
while Patient CLL_125 did not show 11q and 17p mutations. The patient CLL_173 was
observed to have 11q deletion, which is a cytogenetic abnormality. The patients CLL_173
and CLL_175 had relapsing cancers and were previously treated with Fludarabine,
cyclophosphamide and rituximab (FCR) which is a chemotherapeutic drug and a
combination of ofatumumab (OFA)-chlorambucil (CHL)(OFA-Chl) respectively. 

<img width="454" alt="image" src="https://github.com/vanillaberryparfait/University-Of-Glasgow-MSc-Thesis/assets/80147829/71390203-4ebf-4a93-a93c-6d158cf4399d">

These above patients' samples were subjected to bulk RNA-seq analysis under five
distinct growth conditions, which were as follows: <br/>
Group 1 = NTL (non-proliferating) <br/>
Group 2 = NTL-CD40L/IL4 (proliferating) <br/>
Group 3 = NTL-CD40L/IL4 treated with AZD8055 <br/>
Group 4 = NTL-CD40L/IL4 treated with Ibrutinib <br/>
Group 5 = NTL-CD40L/IL4 treated with combination (AZD8055 + IB) <br/>

## Installation

The entire analysis was carried out on RStudio Version 2023.06.0+421 (2023.06.0+421). <br/> 
To run the scripts, the following packages were installed:

```r
install.packages(c("ggplot2", "ggthemes", "ggrepel", "reshape2", "vctrs", "amap", "devtools", "Biobase", "sva", "pheatmap"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DESeq2"))

```

## Workflow

Loading Libraries :
```r
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(reshape2)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(vctrs)
library(amap)
library(DESeq2)
library(devtools)
library(Biobase)
library(sva)
library(pheatmap)
```


Libraries Overview

* ggplot2: A powerful plotting system in R, widely used for creating publication-quality graphics.
* ggthemes: Provides additional themes and options for customizing plots created with ggplot2.
* ggrepel: Ensures labels in ggplot2 plots do not overlap, improving readability.
* reshape2: Useful for transforming and restructuring data frames, facilitating data manipulation.
* clusterProfiler: Performs gene set enrichment analysis, identifying biological pathways enriched in a gene list.
* org.Hs.eg.db: Offers annotation data specific to Homo sapiens, aiding in gene ID conversion and annotation.
* vctrs: Implements methods for creating and manipulating custom vector-like objects.
* amap: Implements clustering and spatial data analysis, useful for geographical and spatial data visualization.
* DESeq2: Performs differential gene expression analysis, identifying genes that are significantly differentially expressed between conditions.
* devtools: Provides tools for R package development, including functions for package installation, documentation, and testing.
* Biobase: Offers fundamental tools and classes for bioinformatics data analysis and manipulation.
* sva: Performs surrogate variable analysis to identify and adjust for hidden sources of variation in high-dimensional data , thereby improving the accuracy of statistical analyses.
* pheatmap: Generates highly customizable heatmaps for visualizing large datasets and patterns within data. 


Loading Gene Count Data and Sample Sets 

```r
g_count_table = read.csv("gene_count.csv",row.names = 1)
ss=read.csv("ss.csv",header=TRUE,row.names = 1, sep=',')
```

Correcting Batch Effect and Differential Expression Analysis
```r
g_count_corrected = ComBat_seq(as.matrix(g_count), batch=ss$PATIENT , group = ss$GROWTHCONDITION)
g_dds=DESeqDataSetFromMatrix(countData=g_count_corrected,colData=ss,design=~ GROWTHCONDITION )
#checks whether g_dds  is a matrix
class(counts(g_dds))
#removes rows less than 1
notAllZero <- (rowSums(counts(g_dds))>0)
g_dds <- g_dds[notAllZero,]
# differential expression (DE) analysis for g_dds 
g_dds= DESeq(g_dds)

```
* Batch Effect:
 The eatch effect was corrected using the ComBat_seq function. Correcting batch effects is essential for improving data quality and accurately interpreting biological insights from  analysis.

* Differential Expression Analysis :
  Differential expression analysis was carried out using the package DESeq2. Prior to
  starting with the differential expression analysis, DESeqDataSetFromMatrix function is
  used to create a DESeq DataSet object. Next, the class function was used to check if the
  DESeq DataSet has been correctly formatted and a matrix has been generated. Rows in
  the matrix having expression less than 1 were further removed as DESeq2 cannot perform
  analysis on genes with such values. Finally, differential expression analysis was carried
  out using the DESeq function. This function performs various steps such as estimating
  size factors, estimating dispersions, calculating gene-wise dispersion estimates,
  examining the mean-dispersion relationship, obtaining final dispersion estimates, and
  fitting a model for testing differential expression. 

Functions Used for the Plots :

```r
## PCA 
em_scaled = na.omit(data.frame(cale(t(em))))
xx = prcomp(t(em_scaled))
pca_coordinates = data.frame(xx$x)

ggp1 = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss$GROWTHCONDITION)) + 
  geom_point(size=6)+ geom_text_repel(aes(label=row.names(ss) ))
png("PCA_GROWTHCONDITION.png", width=1250, height = 1000)
print(ggp1)
ggp1

## Volcano Plots

plot_volcano = function(de_table, p_threshold, fold_threshold ,title) 
{
  #Sort by p value
  sorted = order(de_table[ ,"p"] , decreasing=FALSE) 
  #The following code creates two tables that contain top 10 up regulated and down regulated genes
  de_table = de_table[sorted,]
  de_table_sig = subset(de_table, padj < p_threshold & abs(log2fold) >fold_threshold)
  de_table_sig_up = subset(de_table, padj< p_threshold & log2fold > fold_threshold)
  de_table_sig_down = subset(de_table, padj< p_threshold & log2fold < fold_threshold) 
  de_table_sig_up_top10 = de_table_sig_up[1:10,] 
  de_table_sig_down_top10 = de_table_sig_down[1:10,]
  
  ggp=ggplot(de_table, aes(x = log2fold, y = mlog10p)) + 
    geom_point(colour = "black",shape=15) + 
    geom_point(data = de_table_sig_down, colour = "blue",shape=16) + 
    geom_point(data = de_table_sig_up, colour = "red",shape=18) +
    labs(title=title, x="log2 Fold Change", y="-log10 p-value") + 
    my_theme+ 
    geom_vline(xintercept=-1, linetype="dashed", color = "grey", size=0.5) + 
    geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5) + 
    geom_hline(yintercept=-log10(0.05))+ xlim(c(-20, 20)) + ylim(c(0, 50))+ 
    geom_label_repel(data=de_table_sig_up_top10, aes(label=gene_name), show.legend = TRUE) +
    geom_label_repel(data=de_table_sig_down_top10, aes(label= gene_name),show.legend = TRUE)
  file_name <- paste0(title, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(ggp)
  dev.off() 

  return(ggp) }


## MA plots

plot_MA = function(de_table, p_threshold, fold_threshold,title) 
{
  de_table = na.omit(de_table)
  
  sorted = order(de_table[ ,"p"] , decreasing=FALSE)
  
  de_table = de_table[sorted,]
  de_table_sig = subset(de_table, padj < p_threshold & abs(log2fold) >fold_threshold)
  de_table_sig_up = subset(de_table, padj< p_threshold & log2fold > fold_threshold)
  de_table_sig_down = subset(de_table, padj< p_threshold & log2fold < fold_threshold) 
  de_table_sig_up_top10 = de_table_sig_up[1:10,] 
  de_table_sig_down_top10 = de_table_sig_down[1:10,]
  
  ggp=ggplot(de_table, aes(x = log10(rowMeans), y = log2fold)) + 
    geom_point(colour = "black") + 
    geom_point(data = de_table_sig_down, colour = "blue") + 
    geom_point(data = de_table_sig_up, colour = "red") +
    labs(title=title, x="log10 Mean Expression", y="log2 Fold Change") +
    my_theme + ylim(c(-10,10))+ geom_label_repel(data=de_table_sig_up_top10, aes(label=gene_name)) +
    geom_label_repel(data=de_table_sig_down_top10, aes(label=gene_name))
  file_name <- paste0(title, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(ggp)
  dev.off() 
  
  return(ggp) }

```

* PCA Plot : Essential to visualise clusters and patterns based on similarity from our data and hlp us in understanding data structure and identifying outliers or clusters.
* MA plots : MA plots depict the relationship between the mean expression level (M) and the average abundance (A) of genes or features, highlighting differential expression between    conditions
* Volcano Plots : Volcano Plot is plot between Fold Change and P-value Volcano which helps in indentifying  significantly differentially expressed entities.


Pathway Analysis :

```r
####UP REGULATED PATHWAY RESULTS#### USING GO.
up_pathway_analysis = function(DE){
  
  sig_genes = row.names(subset(DE,padj < 0.05 & abs(log2fold) > 1))
  sig_genes_entrez=bitr(sig_genes, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  pathway_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  return(pathway_results)
}
```

* Over Representation Analysis was carried out  using the
  packages ‘clusterProfiler’ and ‘org.Hs.eg.db’. The function enrichGO was used to carry
  out GO enrichment analysis to understand which biological processes are enriched among
  the differentially expressed genes.
  
Function used to plot a Bar Plot to visualise the Pathway Analysis :

```r

pathway_bar <- function(a, b) {
  barplot <- barplot(a, showCategory = 10)
  file_name <- paste0(b, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(barplot)
  dev.off()  # Close the graphics device
  return(barplot)
}
```

## Overview

Further studies carried out in this study included a more detailed pathway analysis using Metaboanalyst and functional enrichment using DAVID, followed by network analysis with STRING to explore gene correlations and functional enrichments in drug synergy.

Overall, this study investigated the combined effects of AZD8055 and ibrutinib in CLL, uncovering pathways such as AGE-RAGE and FoxO signaling. It revealed insights into cancer biology and potential links to neurodegenerative diseases.







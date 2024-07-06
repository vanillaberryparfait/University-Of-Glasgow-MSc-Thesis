# MSc Thesis Project as Carried Out At The University Of Glasgow 2022-2023

This repository contains the scripts and data used in my MSc thesis project. In this project , I investigated the efficacy of combining the BTK inhibitor ibrutinib and the mTOR inhibitor AZD8055 in treating chronic lymphocytic leukemia (CLL). By analyzing differentially expressed genes in CLL patients under various conditions, I explored synergistic effects, enriched pathways, and potential links between CLL and neurodegenerative disorders

## Table of Contents

- [Data](#data)
- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)

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

## Usage

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



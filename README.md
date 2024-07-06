# MSc Thesis Project as Carried Out At The University Of Glasgow 2022-2023

This repository contains the scripts and data used in my MSc thesis project. In this project , I investigated the efficacy of combining the BTK inhibitor ibrutinib and the mTOR inhibitor AZD8055 in treating chronic lymphocytic leukemia (CLL). By analyzing differentially expressed genes in CLL patients under various conditions, I explored synergistic effects, enriched pathways, and potential links between CLL and neurodegenerative disorders

## Table of Contents

- [Installation](#installation)
- [Data](#data)
- [Usage](#usage)
- [Outputs](#outputs)

## Installation

To run the scripts, you need to have R installed along with the following packages:

```r
install.packages(c("ggplot2", "ggthemes", "ggrepel", "reshape2", "vctrs", "amap", "devtools", "Biobase", "sva", "pheatmap"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DESeq2"))

####LOADING LIBRARIES####
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


####LOADING TABLES#####
#Loading the gene count matrix
g_count_table = read.csv("gene_count.csv",row.names = 1)
ss=read.csv("ss.csv",header=TRUE,row.names = 1, sep=',')
g_count= g_count_table[,c(1:25)]
g_count_table$rowMeans = rowMeans(g_count_table[,1:25])
annotated = g_count_table[,c(26:35)]

ss$GROWTHCONDITION <- factor(ss$GROWTHCONDITION)

#CORRECTING THE BATCH EFFECT
g_count_corrected = ComBat_seq(as.matrix(g_count), batch=ss$PATIENT , group = ss$GROWTHCONDITION)




g_dds=DESeqDataSetFromMatrix(countData=g_count_corrected,colData=ss,design=~ GROWTHCONDITION )
g_dds

#checks whether g_dds  is a matrix
class(counts(g_dds))

#removes rows less than 1
notAllZero <- (rowSums(counts(g_dds))>0)
g_dds <- g_dds[notAllZero,]

#checks how many genes left
g_dds
dim(g_dds)

# differential expression (DE) analysis for g_dds 
g_dds= DESeq(g_dds)

#get em table , corrected
em = round(as.data.frame(counts(g_dds, normalized=TRUE)),2)


## PCA (corrected)
em_scaled = na.omit(data.frame(cale(t(em))))
xx = prcomp(t(em_scaled))
pca_coordinates = data.frame(xx$x)

ggp1 = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss$GROWTHCONDITION)) + 
  geom_point(size=6)+ geom_text_repel(aes(label=row.names(ss) ))
png("PCA_GROWTHCONDITION.png", width=1250, height = 1000)
print(ggp1)
ggp1




### CREATING DE TABLES######


#1vs2
res= results(g_dds,c("GROWTHCONDITION","1","2"))
de =res[order(res$padj),]
DE_1vs2 = as.data.frame(de)
DE_1vs2$ID = row.names(DE_1vs2)
DE_1vs2$mlog10p = -log10(DE_1vs2$pvalue)
DE_1vs2=DE_1vs2[,c(7,2,5,6,8)]
colnames(DE_1vs2)=c("ID","log2fold","p","padj","mlog10p")
DE_1vs2=na.omit(DE_1vs2)
sig_1vs2=subset(DE_1vs2,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_1vs2, file="DE_1vs2.csv",row.names=FALSE, sep="\t", quote = FALSE)

#2vs3
res= results(g_dds,c("GROWTHCONDITION","3","2"))
de =res[order(res$padj),]
DE_3vs2 = as.data.frame(de)
DE_3vs2$ID = row.names(DE_3vs2)
DE_3vs2$mlog10p = -log10(DE_3vs2$pvalue)
DE_3vs2=DE_3vs2[,c(7,2,5,6,8)]
colnames(DE_3vs2)=c("ID","log2fold","p","padj","mlog10p")
DE_3vs2=na.omit(DE_3vs2)
sig_3vs2=subset(DE_3vs2,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_3vs2, file="DE_3vs2.csv",row.names=FALSE, sep="\t", quote = FALSE)

#2vs4
res= results(g_dds,c("GROWTHCONDITION","4","2"))
de =res[order(res$padj),]
DE_2vs4 = as.data.frame(de)
DE_2vs4$ID = row.names(DE_2vs4)
DE_2vs4$mlog10p = -log10(DE_2vs4$pvalue)
DE_2vs4=DE_2vs4[,c(7,2,5,6,8)]
colnames(DE_2vs4)=c("ID","log2fold","p","padj","mlog10p")
DE_2vs4=na.omit(DE_2vs4)
sig_2vs4=subset(DE_2vs4,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_2vs4, file="DE_2vs4.csv",row.names=FALSE, sep="\t", quote = FALSE)

#2vs5
res= results(g_dds,c("GROWTHCONDITION","5","2"))
de =res[order(res$padj),]
DE_2vs5 = as.data.frame(de)
DE_2vs5$ID = row.names(DE_2vs5)
DE_2vs5$mlog10p = -log10(DE_2vs5$pvalue)
DE_2vs5=DE_2vs5[,c(7,2,5,6,8)]
colnames(DE_2vs5)=c("ID","log2fold","p","padj","mlog10p")
DE_2vs5=na.omit(DE_2vs5)
sig_2vs5=subset(DE_2vs5,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_2vs5, file="DE_2vs5.csv",row.names=FALSE, sep="\t", quote = FALSE)

#3vs5
res= results(g_dds,c("GROWTHCONDITION","5","3"))
de =res[order(res$padj),]
DE_3vs5 = as.data.frame(de)
DE_3vs5$ID = row.names(DE_3vs5)
DE_3vs5$mlog10p = -log10(DE_3vs5$pvalue)
DE_3vs5=DE_3vs5[,c(7,2,5,6,8)]
colnames(DE_3vs5)=c("ID","log2fold","p","padj","mlog10p")
DE_3vs5=na.omit(DE_3vs5)
sig_3vs5=subset(DE_3vs5,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_3vs5, file="DE_3vs5.csv",row.names=FALSE, sep="\t", quote = FALSE)

#4vs5
res= results(g_dds,c("GROWTHCONDITION","5","4"))
de =res[order(res$padj),]
DE_4vs5 = as.data.frame(de)
DE_4vs5$ID = row.names(DE_4vs5)
DE_4vs5$mlog10p = -log10(DE_4vs5$pvalue)
DE_4vs5=DE_4vs5[,c(7,2,5,6,8)]
colnames(DE_4vs5)=c("ID","log2fold","p","padj","mlog10p")
DE_4vs5=na.omit(DE_4vs5)
sig_4vs5=subset(DE_4vs5,padj < 0.05 & abs(log2fold) > 1)
write.table(DE_4vs5, file="DE_4vs5.csv",row.names=FALSE, sep="\t", quote = FALSE)


#######creating a dataset for volcano and ma plots #########
de1_2_annotated=merge(DE_1vs2,annotated,by.x=0,by.y=0)
de2_3_annotated=merge(DE_3vs2,annotated,by.x=0,by.y=0)
de2_4_annotated=merge(DE_2vs4,annotated,by.x=0,by.y=0)
de2_5_annotated=merge(DE_2vs5,annotated,by.x=0,by.y=0)
de3_5_annotated=merge(DE_3vs5,annotated,by.x=0,by.y=0)
de4_5_annotated=merge(DE_4vs5,annotated,by.x=0,by.y=0)

####CREATING A THEME####
library(ggthemes)
my_theme <- theme_light() + 
  theme(
    plot.background = element_rect(fill = "lightpink"), 
    plot.title = element_text(size = 15), 
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_line(colour = "gray"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "gray"),
    axis.ticks = element_line(colour = "gray"),
    axis.text = element_text(color = "black"),
    plot.margin=unit(c(1,1,1,1),"cm")
  )
###VOLCANO PLOT#### 
####CREATING VOLCANO PLOTS###
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
###VOLCANO OUTPUT###
volcano_1_2= plot_volcano(de1_2_annotated,0.05, 1,"Volcano plot group 1 vs 2")
volcano_2_3= plot_volcano(de2_3_annotated,0.05, 1,"Volcano plot group 2 vs 3")
volcano_2_4= plot_volcano(de2_4_annotated,0.05, 1,"Volcano plot group 2 vs 4")
volcano_2_5= plot_volcano(de2_5_annotated,0.05, 1,"Volcano plot group 2 vs 5")
volcano_3_5= plot_volcano(de3_5_annotated,0.05, 1,"Volcano plot group 3 vs 5")
volcano_4_5= plot_volcano(de4_5_annotated,0.05, 1,"Volcano plot group 4 vs 5")

### MA PLOT ####
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

##MA OUTPUT####
MA_1_2= plot_MA(de1_2_annotated,0.05, 1,"MA plot group 1 vs 2")
MA_2_3= plot_MA(de2_3_annotated,0.05, 1,"MA plot group 2 vs 3")
MA_2_4= plot_MA(de2_4_annotated,0.05, 1,"MA plot group 2 vs 4")
MA_2_5= plot_MA(de2_5_annotated,0.05, 1,"MA plot group 2 vs 5")
MA_3_5= plot_MA(de3_5_annotated,0.05, 1,"MA plot group 3 vs 5")
MA_4_5= plot_MA(de4_5_annotated,0.05, 1,"MA plot group 4 vs 5")



####Pathway analysis ####

####UP REGULATED PATHWAY RESULTS#### USING GO.
up_pathway_analysis = function(DE){
  
  sig_genes = row.names(subset(DE,padj < 0.05 & abs(log2fold) > 1))
  sig_genes_entrez=bitr(sig_genes, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  pathway_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  return(pathway_results)
}

####BAR PLOT####
pathway_bar <- function(a, b) {
  barplot <- barplot(a, showCategory = 10)
  file_name <- paste0(b, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(barplot)
  dev.off()  # Close the graphics device
  return(barplot)
}

####GO PLOT####
pathway_goplot=function(a,b){
  goplot= goplot(a,showCategory=10)
  file_name <- paste0(b, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(goplot)
  dev.off()  # Close the graphics device
  return(goplot)
}
####DOT PLOT####
pathway_dotplot=function(a,b){
  dotplot= dotplot(a,showCategory=10)
  file_name <- paste0(b, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(dotplot)
  dev.off()  # Close the graphics device
  return(dotplot)
}
#####CNET PLOT####
pathway_cnetplot=function(a,b){
  cnetplot = cnetplot(a, categorySize="p")
  file_name <- paste0(b, ".png")  # Create dynamic file name
  png(file_name, width = 1250, height = 1000)
  print(cnetplot)
  dev.off()  # Close the graphics device
  return(cnetplot)
}

###OUTPUTS: PATHWAY ANALYSIS ###

#1vs2 
up_pathway_analysis_1_2= up_pathway_analysis(DE_1vs2)
bar_up_1_2 = pathway_bar(up_pathway_analysis_1_2,"bar_up_1_2")
go_up_1_2=pathway_goplot(up_pathway_analysis_1_2,"go_up_1_2")
dot_up_1_2=pathway_dotplot(up_pathway_analysis_1_2,"dot_up_1_2")
cnet_up_1_2=pathway_cnetplot(up_pathway_analysis_1_2,"cnet_up_1_2")

#2vs3 
up_pathway_analysis_2_3= up_pathway_analysis(DE_3vs2)
bar_up_2_3 = pathway_bar(up_pathway_analysis_2_3,"bar_up_2_3")
go_up_2_3=pathway_goplot(up_pathway_analysis_2_3,"go_up_2_3")
dot_up_2_3=pathway_dotplot(up_pathway_analysis_2_3,"dot_up_2_3")
cnet_up_2_3=pathway_cnetplot(up_pathway_analysis_2_3,"cnet_up_2_3")


#2vs4
up_pathway_analysis_2_4= up_pathway_analysis(DE_2vs4)
bar_up_2_4 = pathway_bar(up_pathway_analysis_2_4,"bar_up_2_4")
go_up_2_4=pathway_goplot(up_pathway_analysis_2_4,"go_up_2_4")
dot_up_2_4=pathway_dotplot(up_pathway_analysis_2_4,"dot_up_2_4")
cnet_up_2_4=pathway_cnetplot(up_pathway_analysis_2_4,"cnet_up_2_4")

#2vs5
up_pathway_analysis_2_5= up_pathway_analysis(DE_2vs5)
bar_up_2_5 = pathway_bar(up_pathway_analysis_2_5,"bar_up_2_5")
go_up_2_5=pathway_goplot(up_pathway_analysis_2_5,"go_up_2_5")
dot_up_2_5=pathway_dotplot(up_pathway_analysis_2_5,"dot_up_2_5")
cnet_up_2_5=pathway_cnetplot(up_pathway_analysis_2_5,"cnet_up_2_5")


#3vs5
up_pathway_analysis_3_5= up_pathway_analysis(DE_3vs5) 
bar_up_3_5 = pathway_bar(up_pathway_analysis_3_5,"bar_up_3_5")
go_up_3_5=pathway_goplot(up_pathway_analysis_3_5,"go_up_3_5")
dot_up_3_5=pathway_dotplot(up_pathway_analysis_3_5,"dot_up_3_5")
cnet_up_3_5=pathway_cnetplot(up_pathway_analysis_3_5,"cnet_up_3_5")

#4vs5
up_pathway_analysis_4_5= up_pathway_analysis(DE_4vs5)
bar_up_4_5 = pathway_bar(up_pathway_analysis_4_5,"bar_up_4_5")
go_up_4_5=pathway_goplot(up_pathway_analysis_4_5,"go_up_4_5")
dot_up_4_5=pathway_dotplot(up_pathway_analysis_4_5,"dot_up_4_5")
cnet_up_4_5=pathway_cnetplot(up_pathway_analysis_4_5,"cnet_up_4_5")

###arranging emscaled table by groups not patients
emscaled2= em_scaled[,c(1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,14,19,24,5,10,15,20,25)]


########SYNERGY TABLE#########
synergy=read.csv("em_synergy_table.csv", sep=',')
emcopy=em
emcopy$rowname=row.names(emcopy)
row.names(synergy)=emcopy$rowname

### CUTOFF SET AT 0.7 for analysis 
CLL_175=row.names(subset(synergy,(Ratio_CLL_175 > 0.7 )))
CLL_125=row.names(subset(synergy,(Ratio_CLL_125 > 0.7 )))
CLL_203=row.names(subset(synergy,(Ratio_CLL_203 > 0.7 )))
CLL_173=row.names(subset(synergy,(Ratio_CLL_173 > 0.7 )))
CLL_186=row.names(subset(synergy,(Ratio_CLL_186 > 0.7 )))


#####INDIVIDUAL TABLES FOR EACH PATIENT
annotated_copy = annotated
annotated_copy$ID=row.names(annotated_copy)

CLL_175_data=as.data.frame(CLL_175)
names(CLL_175_data)="ID"
CLL_175_annotated=merge(CLL_175_data,annotated_copy,by.y="ID")

CLL_125_data=as.data.frame(CLL_125)
names(CLL_125_data)="ID"
CLL_125_annotated=merge(CLL_125_data,annotated_copy,by.y="ID")

CLL_203_data=as.data.frame(CLL_203)
names(CLL_203_data)="ID"
CLL_203_annotated=merge(CLL_203_data,annotated_copy,by.y="ID")

CLL_173_data=as.data.frame(CLL_173)
names(CLL_173_data)="ID"
CLL_173_annotated=merge(CLL_173_data,annotated_copy,by.y="ID")

CLL_186_data=as.data.frame(CLL_186)
names(CLL_186_data)="ID"
CLL_186_annotated=merge(CLL_186_data,annotated_copy,by.y="ID")


#####PATHWAYS FOR EACH PATIENT , SEPARATELY.
CLL_pathway_analysis = function(DE){
  
  sig_genes_entrez=bitr(DE, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  pathway_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  return(pathway_results)
}

### CLL 175

CLL_175_pathway=CLL_pathway_analysis(CLL_175)

pathway_bar(CLL_175_pathway,"CLL_175_pathway_bar")
pathway_dotplot(CLL_175_pathway,"CLL_175_pathway_dotplot")
pathway_goplot(CLL_175_pathway,"CLL_175_pathway_goplot")
pathway_cnetplot(CLL_175_pathway,"CLL_175_pathway_cnet")

### CLL 125

CLL_125_pathway=CLL_pathway_analysis(CLL_125)

pathway_bar(CLL_125_pathway,"CLL_125_pathway_bar")
pathway_dotplot(CLL_125_pathway,"CLL_125_pathway_dotplot")
pathway_goplot(CLL_125_pathway,"CLL_125_pathway_goplot")
pathway_cnetplot(CLL_125_pathway,"CLL_125_pathway_cnet")


### CLL 203

CLL_203_pathway=CLL_pathway_analysis(CLL_203)

pathway_bar(CLL_203_pathway,"CLL_203_pathway_bar")
pathway_dotplot(CLL_203_pathway,"CLL_203_pathway_dotplot")
pathway_goplot(CLL_203_pathway,"CLL_203_pathway_goplot")
pathway_cnetplot(CLL_203_pathway,"CLL_203_pathway_cnet")

### CLL 173

CLL_173_pathway=CLL_pathway_analysis(CLL_173)

pathway_bar(CLL_173_pathway,"CLL_173_pathway_bar")
pathway_dotplot(CLL_173_pathway,"CLL_173_pathway_dotplot")
pathway_goplot(CLL_173_pathway,"CLL_173_pathway_goplot")
pathway_cnetplot(CLL_173_pathway,"CLL_173_pathway_cnet")

### CLL 186

CLL_186_pathway=CLL_pathway_analysis(CLL_186)

pathway_bar(CLL_186_pathway,"CLL_186_pathway_bar")
pathway_dotplot(CLL_186_pathway,"CLL_186_pathway_dotplot")
pathway_goplot(CLL_186_pathway,"CLL_186_pathway_goplot")
pathway_cnetplot(CLL_186_pathway,"CLL_186_pathway_cnet")


#### ALL PATIENTS, COMMON GENES.
allpatients_common=row.names(subset(synergy,(Ratio_CLL_175 > 0.7 )&(Ratio_CLL_125 > 0.7 )&(Ratio_CLL_203 > 0.7 )&(Ratio_CLL_173 > 0.7)&(Ratio_CLL_186 > 0.7 )))
allpatients_common_data= as.data.frame(allpatients_common)
names(allpatients_common_data)="ID"
allpatients_common_annotated=merge(allpatients_common_data,annotated_copy,by.y="ID")
write.table(allpatients_common_annotated, file="allpatients_common_annotated.csv",row.names=FALSE, sep=",", quote = FALSE)

##PATHWAY ANALYSIS FOR COMMON GENES

all_common_pathway=CLL_pathway_analysis(allpatients_common)

pathway_bar(all_common_pathway,"genes_common_in_all_patients_pathway_bar")
pathway_dotplot(all_common_pathway,"genes_common_in_all_patients_pathway_dotplot")
pathway_goplot(all_common_pathway,"genes_common_in_all_patients_pathway_goplot")
pathway_cnetplot(all_common_pathway,"genes_common_in_all_patients_pathway_cnet")


#####GENERATING TABLES FOR FURTHER ANALYSIS
em2= em[,c(1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,14,19,24,5,10,15,20,25)]
allpatients_common_gene= merge(allpatients_common_annotated,em2,by.x=1,by.y=0)
row.names(allpatients_common_gene)= allpatients_common_gene$ID
allpatients_common_gene = allpatients_common_gene[,c(2,12:36)]
write.table(allpatients_common_gene, file="allpatients_common_gene.csv",row.names=FALSE, sep=",", quote = FALSE)
##In this allpatients_common_gene.csv , a row for groups was added in excel.


##HEATMAP
plot_heatmap=function(signature,title,filename){
  scaled_sig = emscaled2[signature,]
  hm.matrix = as.matrix(scaled_sig)
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  colours = c("blue","white","magenta")
  
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + scale_fill_gradientn(colours = colorRampPalette(colours)(5))+ 
    theme(axis.text.y = element_blank(),plot.background = element_rect(fill = "lightyellow"), axis.ticks=element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'))+
    ggtitle(title)+ geom_vline(xintercept = which(signature %in% rownames(hm.matrix_clustered)), 
                               color = "gray70", linetype = "dotted")
  png(filename, width = 1250, height = 1000)
  print(ggp)
  dev.off() 
  return(ggp)
}

##Heatmap!
sig_heatmapList_1vs2 = sig_1vs2[,c(1)]
sig_heatmapList_1vs2_data=as.data.frame(sig_heatmapList_1vs2)
names(sig_heatmapList_1vs2_data)="ID"
sig_heatmapList_2vs3 = sig_3vs2[,c(1)]
sig_heatmapList_2vs3_data=as.data.frame(sig_heatmapList_2vs3)
names(sig_heatmapList_2vs3_data)="ID"
sig_heatmapList_2vs4 = sig_2vs4[,c(1)]
sig_heatmapList_2vs4_data=as.data.frame(sig_heatmapList_2vs4)
names(sig_heatmapList_2vs4_data)="ID"
sig_heatmapList_2vs5 = sig_2vs5[,c(1)]
sig_heatmapList_2vs5_data=as.data.frame(sig_heatmapList_2vs5)
names(sig_heatmapList_2vs5_data)="ID"

sig_heatmapList_3vs5 = sig_3vs5[,c(1)]
sig_heatmapList_4vs5 = sig_4vs5[,c(1)]

# ALL GROUPS common gene
plot_heatmap(allpatients_common,"ALL GROUPS", "allgroups.png")
#SIG heatmaps
plot_heatmap(sig_heatmapList_1vs2,"heatmap_1vs2","heatmap_1vs2.png")
plot_heatmap(sig_heatmapList_2vs3,"heatmap_2vs3","heatmap_2vs3.png")
plot_heatmap(sig_heatmapList_2vs4,"heatmap_2vs4","heatmap_2vs4.png")
plot_heatmap(sig_heatmapList_2vs5,"heatmap_2vs5","heatmap_2vs5.png")
plot_heatmap(sig_heatmapList_3vs5,"heatmap_3vs5","heatmap_3vs5.png")
plot_heatmap(sig_heatmapList_4vs5,"heatmap_4vs5","heatmap_4vs5.png")

#proliferating vs nonprolif , btk , mtor , comb

all_list1 <- c(sig_heatmapList_1vs2, sig_heatmapList_2vs3, sig_heatmapList_2vs4, sig_heatmapList_2vs5)
plot_heatmap(all_list,"heatmap_2vs1,3,4,5","heatmap_2vs1,3,4,5.png")


## all groups
all_list <- c(sig_heatmapList_1vs2, sig_heatmapList_2vs3, sig_heatmapList_2vs4, sig_heatmapList_2vs5,sig_heatmapList_3vs5,sig_heatmapList_4vs5)
plot_heatmap(all_list,"Heatmap of all significant genes present across all contrasts ","heatmap_all.png")

###heatmap of all synergestic drugs
synergy_all_list <- c(CLL_125, CLL_173, CLL_175, CLL_186,CLL_203)
synergy_all_list = allpatients_common
plot_heatmap(synergy_all_list,"Heatmap of all genes showing synergy ","synergyheatmap_all.png")




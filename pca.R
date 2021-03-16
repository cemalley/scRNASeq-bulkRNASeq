library(data.table)
library(DESeq2)
library(ggplot2)
library(ggrepel)

load('/Volumes/ncatssctl/NGS_related/BulkRNA/GTEx/Seungmi_ISB013/GTEx.Seungmi.DDS.RData')

#how to subset
#dds.subset <- dds[,dds@colData@listData$condition %in% c('H9noc_D28','GTEx_Hypothalamus','GTEx_Hippocampus','GTEx_Cortex',"GTEx_Spinal_cord_cervical_c-1")]
dds.subset <- dds

matrixFile <- as.data.frame(counts(dds.subset, normalized=FALSE))
mat <- as.matrix(matrixFile)

lowCountLimit <- 25
minSamplesLimit <- 3
factorName <- "condition"

normCounts <- data.frame(counts(dds.subset, normalized=TRUE))

normCounts[,"goodCount"]<-rowSums(normCounts > lowCountLimit)

redMat <- subset(normCounts, goodCount >= minSamplesLimit)
redMat <- redMat[,1:(ncol(redMat)-1)]
conds <- factor(dds.subset@colData$condition)

pca <- prcomp(t(redMat), center=TRUE, scale=TRUE)
conditions <- data.frame("Condition"= conds)

pca.df <- as.data.frame(pca$x)

ggplot(pca.df, aes(x=pca$x[,"PC1"], y=pca$x[,"PC2"])) + geom_point(aes(color=conditions$Condition), size=5, alpha=0.5) + theme_bw() +  #geom_text_repel(aes(label=names(pca$x[,"PC1"])),color="black") +
  labs(x="PC1", y="PC2", title="PCA of BrainSpan samples plus Seungmi's cerebral organoids") +
  guides(color=guide_legend(title="Condition"))+theme(panel.grid=element_line(size=1), axis.text = element_text(size=12), axis.title=element_text(face='plain'), title=element_text(face = 'bold'), legend.text = element_text(size=12), legend.title = element_text(face='plain', size=12))



#+ scale_color_manual(values=c('GTEx_Amygdala'='#1f78b4', 'GTEx_Anterior_cingulate_cortex_BA24'='#08306b', 'GTEx_Caudate_basal_ganglia'='#dd3497', 'GTEx_Cerebellar_Hemisphere'='#fb9a99', 'GTEx_Cerebellum'='#e31a1c', 'GTEx_Cortex'='#fdbf6f', 'GTEx_Frontal_Cortex_BA9'='#ff7f00', 'GTEx_Hippocampus'='#cab2d6', 'GTEx_Hypothalamus'='#6a3d9a', 'GTEx_Nucleus_accumbens_basal_ganglia'='#ffff99', 'GTEx_Putamen_basal_ganglia'='#b15928', 'GTEx_Spinal_cord_cervical_c-1'='#8dd3c7', 'GTEx_Substantia_nigra'='#525252', 'H9ES_P31_Y_d60'='#addd8e', 'H9ES_P31_CE_d60'='#41ab5d', 'H9ES_P31_CEPT_d60'='#006837'))


# 3D PCA-----
install.packages('pca3d')
install.packages('rgl')
library(pca3d)
data(metabo)
pca <- prcomp(metabo[,-1], scale.=TRUE)
gr <- factor(metabo[,1])
summary(gr)


#

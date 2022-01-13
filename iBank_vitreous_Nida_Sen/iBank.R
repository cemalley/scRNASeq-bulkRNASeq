library(data.table)
library(Seurat)
library(stringr)

setwd("~/scRNASeq/iBank_vitreous/Countfiles_ENSG")

files <- Sys.glob('*csv')
for (file in files){
  dt <- as.data.frame(fread(file))
  names(dt)[1] <- 'ENSG'
  dt <- as.data.table(dt)
  ensg.genes <- data.table("ENSG" = dt$ENSG)
  genes <- ensg.genes$ENSG
  G_list <- fread('~/BulkRNASeq/Biomart_hsapiens_ensembl_gene_symbols.csv')
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG %in% ensg.genes$ENSG)
  dt <- merge(dt, ensg.genes, by="ENSG")
  dt <- dt[,ENSG:=NULL]
  dt <- dt[!duplicated(external_gene_name),]
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  
  sample_id <- str_split_fixed(file, "_dense_expression_matrix.csv",2)[1]
  
  fwrite(dt, paste0("~/scRNASeq/iBank_vitreous/Counfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}



reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- x[external_gene_name != "NA",]
  #  x <- subset(x, select=c("external_gene_name", unlist(names(x))[2:(length(x))] ))
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x) <- barcodes
  x <- CreateSeuratObject(counts = x, project = "iBank_vitreous")
  x@meta.data$sample <- samplename
  return(x)
}

setwd('~/scRNASeq/iBank_vitreous/Counfiles_gene_symbol/')
files <- Sys.glob('*csv')


iBank_137 <- reformat_for_seurat(fread(files[1]), "Vitreous_137")
iBank_304 <- reformat_for_seurat(fread(files[2]), "Vitreous_304")
iBank_312 <- reformat_for_seurat(fread(files[3]), "Vitreous_312")
iBank_321 <- reformat_for_seurat(fread(files[4]), "Vitreous_321")

finish_seurat <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x, features = row.names(as.data.frame(x@assays$RNA@data)))
  Idents(x) <- 'sample'
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x)
  x <- RunTSNE(x) # error for Noci sample.
  x <- RunUMAP(x, reduction='pca', dims=1:20)
  return(x)
}


iBank_137 <- finish_seurat(iBank_137)
iBank_304 <- finish_seurat(iBank_304)
iBank_312 <- finish_seurat(iBank_312)
iBank_321 <- finish_seurat(iBank_321)

# december 16, 2021-----
setwd('~/scRNASeq/iBank_vitreous/')
load('iBank_304.Seurat.RData')

TSNEPlot(iBank_304)
UMAPPlot(iBank_304)

iBank_304 <- FindClusters(iBank_304, resolution = 0.5)
UMAPPlot(iBank_304)

DE_304 <- FindAllMarkers(iBank_304, only.pos = F, return.thresh = 0.001)
DE_304 <- as.data.table(DE_304)
DE_304 <- DE_304[!grep('^MT-', gene),]
DE_304 <- as.data.table(DE_304)
DE_304 <- DE_304[p_val_adj <= 0.001,]
fwrite(DE_304, 'iBank_304_DE.csv')
DE_304_top <- setorder(setDT(DE_304), -avg_log2FC)[, head(.SD, 5), keyby = 'cluster']

DoHeatmap(iBank_304, features=DE_304_top$gene, raster=F) + guides(color='none')+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color='none') + theme(axis.text.y = element_text(size=8))

iBank_304 <- FindClusters(iBank_304, resolution=0.1)
UMAPPlot(iBank_304)


mm3 <- fread('mm3_list.txt', header=F)

DE_304_mm3 <- FindAllMarkers(iBank_304, features = intersect(rownames(iBank_304), mm3$V1), return.thresh = 0.05)
DE_304_mm3 <- as.data.table(DE_304_mm3)
DE_304_mm3<- DE_304_mm3[p_val_adj <= 0.05,]
DE_304_mm3 <- DE_304_mm3[order(-abs(avg_log2FC))]
DE_304_mm3 <- DE_304_mm3[abs(avg_log2FC)>1,]


DoHeatmap(iBank_304, features=c(DE_304_mm3$gene), raster=F, angle=0, hjust=0) +
  guides(color='none')+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )+
  theme(axis.text.y = element_text(size=8))+labs(title='iBank 304 - genes from literature where abs(avg log2fc)>1')

######

DE_304_mm3$gene[DE_304_mm3$gene %in% DE_304[p_val_adj <= 0.05 & abs(avg_log2FC)>=1,gene]]# 54 of 527


litgenes <- as.data.table(read_excel('literature_genes.xlsx'))

DoHeatmap(iBank_304, features=c(litgenes$`RA FCRL4+`), raster=F, angle=0, hjust=0) +
  guides(color='none')+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )+
  theme(axis.text.y = element_text(size=8))+
  labs(title='iBank 304 - genes from literature\nRA FCRL4+ cell genes')





DE_304 <- FindAllMarkers(iBank_304, return.thresh = 0.001)

# qc-------
load('iBank_304.Seurat.RData')
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
iBank_304[["percent.mt"]] <- PercentageFeatureSet(iBank_304, pattern = "^MT-")
Idents(iBank_304) <- 'sample'
# Visualize QC metrics as a violin plot
VlnPlot(iBank_304, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(iBank_304, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(iBank_304, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

iBank_304 <- subset(iBank_304, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

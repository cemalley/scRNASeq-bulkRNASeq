library(Seurat)
library(dplyr)
library(data.table)
aggr6 <- fread("Aggregation6.allLA_day0_day1_day2_day3_day5_day7.csv") #DDSeq sparse data matrix
aggr6 <- aggr6[!grep("^RP", GeneId),] # remove ribosomal proteins

pdata <- fread("Aggregation6.pDataCDS.csv") # phenotype data table from Monocle analysis on the same data

#Monocle "state" is a segment of the pseudotime trajectory
state1.barcodes <- pdata[State=="1",c("Barcode")]
state2.barcodes <- pdata[State=="2",c("Barcode")]
state3.barcodes <- pdata[State=="3",c("Barcode")]
state4.barcodes <- pdata[State=="4",c("Barcode")]
state5.barcodes <- pdata[State=="5",c("Barcode")]
state6.barcodes <- pdata[State=="6",c("Barcode")]
state7.barcodes <- pdata[State=="7",c("Barcode")]

state3.day3 <- pdata[State=="3" & Sample=="day3",c("Barcode")]
state6.day3 <- pdata[State=="6" & Sample=="day3",c("Barcode")]

aggr6.st1 <- subset(aggr6, select=c(names(aggr6) %in% state1.barcodes$Barcode))
aggr6.st2 <- subset(aggr6, select=c(names(aggr6) %in% state2.barcodes$Barcode))
aggr6.st3 <- subset(aggr6, select=c(names(aggr6) %in% state3.barcodes$Barcode))
aggr6.st4 <- subset(aggr6, select=c(names(aggr6) %in% state4.barcodes$Barcode))
aggr6.st5 <- subset(aggr6, select=c(names(aggr6) %in% state5.barcodes$Barcode))
aggr6.st6 <- subset(aggr6, select=c(names(aggr6) %in% state6.barcodes$Barcode))
aggr6.st7 <- subset(aggr6, select=c(names(aggr6) %in% state7.barcodes$Barcode))

aggr6.st3.day3 <- subset(aggr6, select=c(names(aggr6) %in% state3.day3$Barcode))
aggr6.st6.day3 <- subset(aggr6, select=c(names(aggr6) %in% state6.day3$Barcode))


aggr6.st1 <- as.data.table(aggr6.st1)
aggr6.st2 <- as.data.table(aggr6.st2)
aggr6.st3 <- as.data.table(aggr6.st3)
aggr6.st4 <- as.data.table(aggr6.st4)
aggr6.st5 <- as.data.table(aggr6.st5)
aggr6.st6 <- as.data.table(aggr6.st6)
aggr6.st7 <- as.data.table(aggr6.st7)

aggr6.st3.day3 <- as.data.table(aggr6.st3.day3)
aggr6.st6.day3 <- as.data.table(aggr6.st6.day3)

aggr6 <- as.data.frame(aggr6)


# DE function for repeating with many comparisons-----

DErun <- function(x, y, varname1, varname2){
  var1 <- as.data.table(x)
  var1 <- as.data.frame(var1)
  row.names(var1)<- aggr6$GeneId
  var1 <- CreateSeuratObject(raw.data = var1, project = paste0("NPC.", varname1), min.cells = 5, min.genes = 200)
  var1@meta.data$sample <- varname1
  var1 <- FilterCells(var1, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
  var1 <- NormalizeData(var1)
  var1 <- ScaleData(var1, display.progress = F)
  
  var2 <- as.data.table(y)
  var2 <- as.data.frame(var2)
  row.names(var2)<- aggr6$GeneId
  var2 <- CreateSeuratObject(raw.data = var2, project = paste0("NPC.", varname2), min.cells = 5, min.genes = 200)
  var2@meta.data$sample <- "var2"
  var2 <- FilterCells(var2, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
  var2 <- NormalizeData(var2)
  var2 <- ScaleData(var2, display.progress = F)
  
  seurat.obj <- MergeSeurat(var1, var2, project="NPC", min.cells = 5, min.genes = 200)
  seurat.obj <- SetAllIdent(object = seurat.obj, id = "sample")
  
  mito.genes <- grep(pattern = "^MT", x = rownames(x = seurat.obj@data), value = TRUE)
  percent.mito <- Matrix::colSums(seurat.obj@raw.data[mito.genes, ])/Matrix::colSums(seurat.obj@raw.data)
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = percent.mito, col.name = "percent.mito")
  seurat.obj <- FilterCells(object = seurat.obj, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(25000))
  seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  par(mfrow = c(1, 1))
  seurat.obj <- FindVariableGenes(object = seurat.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  seurat.obj <- ScaleData(object = seurat.obj, vars.to.regress = c("nUMI"))
  seurat.obj <- RunPCA(object = seurat.obj, pc.genes = seurat.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  seurat.obj <- JackStraw(object = seurat.obj, num.replicate = 100)
  seurat.obj <- FindClusters(object = seurat.obj, reduction.type = "pca", dims.use = 1:12, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  seurat.obj <- SetAllIdent(object = seurat.obj, id = "sample")
  markers <- FindMarkers(object = seurat.obj, ident.1 = varname1, min.pct = 0.25)
  markers$GeneId <- row.names(markers)
  names(markers) <- c("p_val", "avg_logFC", paste0("percent.",varname1), paste0("percent.", varname2), "p_val_adj", "GeneId")
  return(markers)
}

markers <- DErun(aggr6.st3.day3, aggr6.st6.day3, "state3.day3", "state6.day3")
varname1 <- "state3.day3"
varname2 <- "state6.day3"
fwrite(markers, paste0("DEbystate/", varname1, ".vs.", varname2, ".markers.noRP.csv"), sep="\t", col.names = T, row.names = F)

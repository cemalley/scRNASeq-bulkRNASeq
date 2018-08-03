library(DESeq2)
library(RColorBrewer)
library(gplots)
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt", suppressUpdates = TRUE)
library(biomaRt)
library(data.table)
library(stringr)

# convert ENSG.# genes to ENSG, then gene symbols----

setwd("Countfiles_ENSG") # htseq count file output from Chromium data comes with ENSG gene IDs
files <- Sys.glob("Sample*counts.txt")

for (file in files){
  
  dt <- fread(file)
  dt <- dt[1:(nrow(dt)-5),]
  names(dt) <- c("ENSG.full", "Counts")
  
  ensg.genes <- data.table("ENSG.full" = dt$ENSG.full)
  
  ensg.genes[,ENSG.short := tstrsplit(ENSG.full, "\\.")[1]]
  
  genes <- ensg.genes$ENSG.short
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG.short", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG.full %in% ensg.genes$ENSG.full)
  dt <- merge(dt, ensg.genes, by="ENSG.full")
  dt <- dt[,c("external_gene_name", "Counts")]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_htseq_counts.txt",2)[1]
  
  fwrite(dt, paste0("Countfiles_gene_symbols/", sample_id, "_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
  
}

# start deseq2 work----

sampleTable <- fread("sampletable.csv")
sampleFiles <- gtools::mixedsort(grep("Sample_",list.files("./Countfiles_gene_symbols/"),value=TRUE))
sampleCondition <- unlist(sampleTable[,condition])
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable
directory <- "Countfiles_gene_symbols/"

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)



# DE TESTS----
res <- results(dds)
res.day0.vs.day4 <- lfcShrink(dds, contrast=c("condition","H9noc_D4", "H9noc_D0")) # default padjust method is Benjamini-Hochberg. must specify "later timepoint" versus "early timepoint" to be in the same log fold change direction as "early versus later" in results.
res.day0.vs.day8 <- lfcShrink(dds, contrast=c("condition","H9noc_D8", "H9noc_D0"))
res.day0.vs.day12 <- lfcShrink(dds, contrast=c("condition","H9noc_D12", "H9noc_D0"))
res.day0.vs.day21 <- lfcShrink(dds, contrast=c("condition","H9noc_D21", "H9noc_D0"))
res.day0.vs.day28 <- lfcShrink(dds, contrast=c("condition","H9noc_D28", "H9noc_D0"))

res.day4.vs.day8 <- lfcShrink(dds, contrast=c("condition","H9noc_D8", "H9noc_D4"))
res.day4.vs.day12 <- lfcShrink(dds, contrast=c("condition","H9noc_D12", "H9noc_D4"))
res.day4.vs.day21 <- lfcShrink(dds, contrast=c("condition","H9noc_D21", "H9noc_D4"))
res.day4.vs.day28 <- lfcShrink(dds, contrast=c("condition","H9noc_D28", "H9noc_D4"))

res.day8.vs.day12 <- lfcShrink(dds, contrast=c("condition","H9noc_D12", "H9noc_D8"))
res.day8.vs.day21 <- lfcShrink(dds, contrast=c("condition","H9noc_D21", "H9noc_D8"))
res.day8.vs.day28 <- lfcShrink(dds, contrast=c("condition","H9noc_D28", "H9noc_D8"))

res.day12.vs.day21 <- lfcShrink(dds, contrast=c("condition","H9noc_D21", "H9noc_D12"))
res.day12.vs.day28 <- lfcShrink(dds, contrast=c("condition","H9noc_D28", "H9noc_D12"))

res.day21.vs.day28 <- lfcShrink(dds, contrast=c("condition","H9noc_D28", "H9noc_D21"))

setwd(directory)
save(dds, res.day0.vs.day12, res.day0.vs.day21, res.day0.vs.day28, res.day0.vs.day4, res.day0.vs.day8, res.day12.vs.day21, res.day12.vs.day28, res.day21.vs.day28, res.day4.vs.day12, res.day4.vs.day21, res.day4.vs.day28, res.day4.vs.day8, res.day8.vs.day12, res.day8.vs.day21, res.day8.vs.day28, sampleTable, file="Tao_nociceptor_DDS.Rdata")

#MA plots----
plotMA(res.day0.vs.day28, alpha=(1*10^(-5)))

#reformat to data.table and add mean per level----
reformat.res <- function(x, condition1, condition2){
  x <- as.data.frame(x)
  x$GeneId <- row.names(x)
  x <- as.data.table(x)
  x[,baseMean:=NULL]
  x[,baseMean1 := rowMeans(counts(dds,normalized=TRUE)[,dds$condition == condition1])]
  x[,baseMean2 := rowMeans(counts(dds,normalized=TRUE)[,dds$condition == condition2])]
  x <- na.omit(x)
  x <- x[order(rank(padj))]
  x <- x[,c("GeneId", "baseMean1", "baseMean2", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  names(x)[2:3] <- c(condition1, condition2)
  return(x)
}

res.day0.vs.day4 <- reformat.res(res.day0.vs.day4, "H9noc_D0", "H9noc_D4")
res.day0.vs.day8 <- reformat.res(res.day0.vs.day8, "H9noc_D0", "H9noc_D8")
res.day0.vs.day12 <- reformat.res(res.day0.vs.day12, "H9noc_D0", "H9noc_D12")
res.day0.vs.day21 <- reformat.res(res.day0.vs.day21, "H9noc_D0", "H9noc_D21")
res.day0.vs.day28 <- reformat.res(res.day0.vs.day28, "H9noc_D0", "H9noc_D28")

res.day4.vs.day8 <- reformat.res(res.day4.vs.day8, "H9noc_D4", "H9noc_D8")
res.day4.vs.day12 <- reformat.res(res.day4.vs.day12, "H9noc_D4", "H9noc_D12")
res.day4.vs.day21 <- reformat.res(res.day4.vs.day21, "H9noc_D4", "H9noc_D21")
res.day4.vs.day28 <- reformat.res(res.day4.vs.day28, "H9noc_D4", "H9noc_D28")

res.day8.vs.day12 <- reformat.res(res.day8.vs.day12, "H9noc_D8", "H9noc_D12")
res.day8.vs.day21 <- reformat.res(res.day8.vs.day21, "H9noc_D8", "H9noc_D21")
res.day8.vs.day28 <- reformat.res(res.day8.vs.day28, "H9noc_D8", "H9noc_D28")

res.day12.vs.day21 <- reformat.res(res.day12.vs.day21, "H9noc_D12", "H9noc_D21")
res.day12.vs.day28 <- reformat.res(res.day12.vs.day28, "H9noc_D12", "H9noc_D28")

res.day21.vs.day28 <- reformat.res(res.day21.vs.day28, "H9noc_D21", "H9noc_D28")

# write----
setwd("DE")
fwrite(res.day0.vs.day4, "DE.H9noc_D0_vs_D4.csv") #
fwrite(res.day0.vs.day8, "DE.H9noc_D0_vs_D8.csv")
fwrite(res.day0.vs.day12, "DE.H9noc_D0_vs_D12.csv")
fwrite(res.day0.vs.day21, "DE.H9noc_D0_vs_D21.csv")
fwrite(res.day0.vs.day28, "DE.H9noc_D0_vs_D28.csv")
fwrite(res.day4.vs.day8, "DE.H9noc_D4_vs_D8.csv")#
fwrite(res.day4.vs.day12, "DE.H9noc_D4_vs_D12.csv")
fwrite(res.day4.vs.day21, "DE.H9noc_D4_vs_D21.csv")
fwrite(res.day4.vs.day28, "DE.H9noc_D4_vs_D28.csv")
fwrite(res.day8.vs.day12, "DE.H9noc_D8_vs_D12.csv")#
fwrite(res.day8.vs.day21, "DE.H9noc_D8_vs_D21.csv")
fwrite(res.day8.vs.day28, "DE.H9noc_D8_vs_D28.csv")
fwrite(res.day12.vs.day21, "DE.H9noc_D12_vs_D21.csv")#
fwrite(res.day21.vs.day28, "DE.H9noc_D21_vs_D28.csv")#

#volcano plot----

#restoplot <- res.day0.vs.day4 # used 3 for lfc threshold
#restoplot <- res.day4.vs.day8 # 6
#restoplot <- res.day8.vs.day12 # 5
#restoplot <- res.day0.vs.day28 
#restoplot <- res.day12.vs.day21 #5
restoplot <- res.day21.vs.day28

alpha <- 1*10^(-50) # Threshold on the adjusted p-value
cols <- densCols(restoplot$log2FoldChange, -log10(restoplot$pvalue))
plot(restoplot$log2FoldChange, -log10(restoplot$padj), col=cols, panel.first=grid(),
     main="Volcano plot: D21 vs D28", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6, xlim=c(-10,10))
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(restoplot$log2FoldChange) > 4 & restoplot$padj < alpha 
text(restoplot$log2FoldChange[gn.selected],
     -log10(restoplot$padj)[gn.selected],
     lab=(restoplot$GeneId)[gn.selected ], cex=0.7)

# condition1 higher is on the left, condition2 higher is on the right

# for plotting counts of a particular gene-----
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene="HECW1", intgroup="condition")

# prepare for GO tests----
setwd("GO/")

filtered <- res.day0.vs.day4[abs(log2FoldChange)>1 & padj < 0.001 & H9noc_D0 > 25 & H9noc_D4 > 25,] # filtered set.
filtered <- as.data.table(filtered$GeneId[1:200])
fwrite(filtered, "H9noc_D0_vs_D4.txt", col.names = F, row.names = F)

filtered <- res.day4.vs.day8[abs(log2FoldChange)>1 & padj < 0.001 & H9noc_D4 > 25 & H9noc_D8 > 25,] # filtered set.
filtered <- as.data.table(filtered$GeneId[1:200])
fwrite(filtered, "H9noc_D4_vs_D8.txt", col.names = F, row.names = F)

filtered <- res.day8.vs.day12[abs(log2FoldChange)>1 & padj < 0.001 & H9noc_D8 > 25 & H9noc_D12 > 25,] # filtered set.
filtered <- as.data.table(filtered$GeneId[1:200])
fwrite(filtered, "H9noc_D8_vs_D12.txt", col.names = F, row.names = F)

filtered <- res.day12.vs.day21[abs(log2FoldChange)>1 & padj < 0.001 & H9noc_D12 > 25 & H9noc_D21 > 25,] # filtered set.
filtered <- as.data.table(filtered$GeneId[1:200])
fwrite(filtered, "H9noc_D12_vs_D21.txt", col.names = F, row.names = F)

filtered <- res.day21.vs.day28[abs(log2FoldChange)>1 & padj < 0.001 & H9noc_D21 > 25 & H9noc_D28 > 25,] # filtered set.
filtered <- as.data.table(filtered$GeneId[1:200])
fwrite(filtered, "H9noc_D21_vs_D28.txt", col.names = F, row.names = F)




#

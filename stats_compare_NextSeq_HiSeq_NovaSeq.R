library(data.table)


raslseq.stats <- as.data.table(read_xlsx('/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/raslseq10_counts.xlsx'))
raslseq.stats #1629 genes. yield: 28.34 Gbp. cycles: 50 | 8 | 8.

load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/ISB015.DDS.RData")

ISB015 <- as.data.frame(as.matrix(counts(dds, normalized=F)))
ISB015.colsums <- colSums(ISB015) #34292 genes
mean(ISB015.colsums)

load("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.Seurat.RData")
#IS019 in house, nextseq

IS019 <- as.data.frame(as.matrix(is019@assays$RNA@counts))
IS019.colsums <- colSums(IS019)
mean(IS019.colsums) # 4342.317 per cell

# is022 in house, novaseq


library(data.table)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(gtools)
library(DESeq2)
library(scales)

load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/VJ4001_D0-30_astrocyte_DDS.RData")

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

genes <- c('SALL4', 'DUSP6', 'NANOG', 'ESRG', 'PAX6', 'SIX3', 'SIX6', 'HESX1', 'FABP7', 'GLAST', 'TNC', 'VIM', 'S100B', 'GLUL', 'SOX9', 'NFIA', 'NFIB', 'CD44', 'GFAP', 'GJA1', 'AQP4', 'ALDH1L1', 'NEUROG1', 'NEUROG2', 'TUBB3', 'SYN1', 'PLP1', 'MOG', 'MBP', 'CLDN5', 'ITM2A', 'ESAM', 'CX3CR1', 'CCL3', 'TNF', 'MKI67')

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genes]
length(found)
length(genes)
genes[!genes %in% found] #[1] "GLAST"

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)

names(mat)<- gsub('_',' ',coldata$condition)

mat

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

h1 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

h1
# 
# myrowanno <- as.data.frame(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Genes for the temporal bulkRNA seq heatmap.xlsx', sheet=2))
# row.names(myrowanno) <- myrowanno$GeneId
# myrowanno <- subset(myrowanno, select=c('Geneset'))
# myrowanno
# 
# ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'), col=list(Geneset=c('Pluripotency'='#a6cee3', 'Neuroepithelium'='#1f78b4', 'Radial Glial'='#b2df8a', 'Astrocyte Precursors'='#33a02c', 'Astrocytes'='#fb9a99', 'Neurons'='#e31a1c', 'Oligodendrocytes'='#fdbf6f', 'Endothelial'='#ff7f00', 'Microglia'='#cab2d6', 'Proliferation'='#6a3d9a')))
# 
# #draw(ha, 1:37)
# 
# draw(h1 + ha, column_title = "Astrocytes D0 - D30")
# 
# 
# h2 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
#               column_names_side = "top",col = cols.use, show_column_names = T,
#               cluster_rows = FALSE, cluster_columns = FALSE,
#               heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))
# 
# draw(h2 + ha, column_title = "Astrocytes D0 - D30")

# rearranged genes to match the genelist
# mat <- read.csv('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/heatmap_mat.csv', row.names = 1, as.is = T, check.names = F)
# 
# h3 <- Heatmap(as.matrix(mat), row_names_side = "right",
#               column_names_side = "top",col = cols.use, show_column_names = T,
#               cluster_rows = FALSE, cluster_columns = FALSE,
#               heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))
# draw(h3+ha, column_title='Astrocytes D0 - D30')
# 
# h4 <- Heatmap(as.matrix(mat), row_names_side = "right",
#               column_names_side = "top",col = cols.use, show_column_names = T,
#               cluster_rows = FALSE, cluster_columns = FALSE,
#               heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))
# 
# draw(h4+ha, column_title='Astrocytes D0 - D30')

library(monocle)
library(dplyr)
library(data.table)

# aggregation 6 - fixed LA_day2 and LA_day5 sample labelling switch-----
exp <- as.data.frame(fread("Aggregation6.allLA_day0_day1_day2_day3_day5_day7.csv")) # genes x counts
exp.samples <- names(exp)[-1]
ss <- data.frame("Barcode"= exp.samples, "Sample"= c(rep("day0", 198), rep("day1", 280), rep("day2", 214), rep("day3", 212), rep("day5", 195), rep("day7", 286))) 

ann <- data.frame("Gene"=exp[1], "gene_short_name"=exp[1]) #gene names - only have hgnc_symbol
exp <- exp[,2:(dim(exp)[2])]
names(ann)<- c("Gene", "gene_short_name")
rownames(ss)<-ss$Barcode

# monocle run--------

pd <- new("AnnotatedDataFrame", data = ss)
fd <- new("AnnotatedDataFrame", data = ann)

cds <- newCellDataSet(as.matrix(exp), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 50))

# QC steps ----
pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))

cds <- cds[,pData(cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) +
                     2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) -
                     2*sd(log10(pData(cds)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(cds), color = Sample, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
             pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)

# Log-transform each value in the expression matrix.
L <- log(exprs(cds[expressed_genes,]))

# clustering cells without marker genes----
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 9, reduction_method = 'tSNE', verbose = T)

cds <- clusterCells(cds, num_clusters = 5) # must change num_clusters for each aggregation run of monocle
plot_cell_clusters(cds, color = "Sample")

# order genes in pseudotime

diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~Sample")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)


#pseudotime trajectory building----
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Sample")

cds_expressed_genes <-  row.names(subset(fData(cds), num_cells_expressed >= 10))

cds_filtered <- cds[cds_expressed_genes,]

# differential gene expression test----
clustering_DEG_genes <- differentialGeneTest(cds_filtered, fullModelFormulaStr = '~Sample', cores = 3) # differential gene expression tested with Sample in the model
clustering_DEG_genes %>% arrange(qval) %>% head()
fwrite((clustering_DEG_genes %>% arrange(qval)), "Aggregation6.clusteringDEGgenes.csv", sep="\t", col.names = T, row.names = F, quote=F)

#trajectories-----
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "Hours")
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 2)
plot_cell_trajectory(cds, color_by = "num_genes_expressed")
plot_cell_trajectory(cds, color_by = "Total_mRNAs")

# coloring by gene expression or any metadata column ------
aggr6 <- fread("Aggregation6.allLA_day0_day1_day2_day3_day5_day7_SWITCH-FIXEDApr3.csv")
pou5f1.data <- aggr6[GeneId=="POU5F1",c(2:1386)]
cells.keep <- row.names(pData(cds))
pou5f1.data <- subset(pou5f1.data, select=c(names(pou5f1.data) %in% cells.keep))
pou5f1.data <- unlist(c(pou5f1.data), use.names = F)
pou5f1.data <- log2(pou5f1.data +1)

pData(cds)$POU5F1 <- pou5f1.data

trajectory.plot <- plot_cell_trajectory(cds, color_by = "POU5F1")
trajectory.plot + scale_color_gradient(low = "#ffcccc", high="red")

# color by gene expression without having to add metadata---
plot_cell_clusters(cds, markers=c("PAX6", "SOX11", "POU5F1"), cell_size = 0.5) + scale_color_gradient(low="green", high="red")

plot_cell_clusters(cds, markers=c("PAX6"), cell_size = 1) + scale_color_gradient(low="green", high="red")
plot_cell_clusters(cds, markers=c("SOX11"), cell_size = 1) + scale_color_gradient(low="green", high="red")
plot_cell_clusters(cds, markers=c("POU5F1"), cell_size = 1) + scale_color_gradient(low="green", high="red")

markerlist <- c("NANOG", "FZD5", "FOXG1", "HES4", "MT2A", "FOXG1", "PRSS23", "GADD45A", "DLK1")
plot_cell_clusters(cds, markers=c(markerlist), cell_size = 1) + scale_color_gradient(low="green", high="red")

# coloring by cell cycle phase.-----
cell.cycle <- fread("Cell.cycle.by.pseudotime.csv")

cells.keep <- row.names(pData(cds))
cells.keep <- as.data.table(cells.keep)
cells.keep$order <- 1:1320

cell.cycle <- merge(cell.cycle, cells.keep, by.x="barcode", by.y="cells.keep")
cell.cycle<-cell.cycle[order(order)] 
cell.cycle.vector <- unlist(c(cell.cycle$phase), use.names=F)

pData(cds)$Cell.cycle <- cell.cycle.vector

plot_cell_trajectory(cds, color_by = "Cell.cycle", cell_size = 0.7) + scale_color_manual(values = c("pink", "red", "orange", "lightgrey", "skyblue"))

fwrite(pData(cds), "Aggregation6.pDataCDS.csv", sep="\t", col.names = T, row.names=T, quote=F)
head(pData(cds_filtered)) #pseudotime is now a column in the dataset.

# states. each state is simply a section of the pseudotime tree between nodes-----
state1 <- row.names(subset(pData(cds), pData(cds)$State=="1"))
state2 <- row.names(subset(pData(cds), pData(cds)$State=="2"))
state3 <- row.names(subset(pData(cds), pData(cds)$State =="3"))
state4 <- row.names(subset(pData(cds), pData(cds)$State =="4"))
state5 <- row.names(subset(pData(cds), pData(cds)$State =="5"))
state6<- row.names(subset(pData(cds), pData(cds)$State=="6"))
state7<- row.names(subset(pData(cds), pData(cds)$State=="7"))

state1 <- as.data.table(state1)
state2 <- as.data.table(state2)
state3 <- as.data.table(state3)
state4 <- as.data.table(state4)
state5 <- as.data.table(state5)
state6 <- as.data.table(state6)
state7 <- as.data.table(state7)


monocle.data.pseudotime <- as.data.frame(pData(cds_filtered))
monocle.data.pseudotime$barcode <- row.names(monocle.data.pseudotime)
monocle.data.pseudotime <- as.data.table(monocle.data.pseudotime)
#Saved Aggregation6-monocle.RData
#load("Aggregation6-monocle.RData")

my_pseudotime_de <- differentialGeneTest(cds_filtered, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 3) #run pseudotime DE test to find genes that change along pseudotime and a function of known Sample (day)

my_pseudotime_de %>% arrange(qval) %>% head()

#plot in pseudotime with markers ----
my_pseudotime_gene <- c("NANOG", "POU5F1", "SOX2", "MAP2", "OTX2", "PAX6", "GBX2")
my_pseudotime_gene <- c("CCNA2", "CCNB1", "CCNC", "CCND2", "CCND3", "CCNE1", "CCNF", "CCNG1", "CCNG2", "CCNH", "CCNT1", "CCNT2")
my_pseudotime_gene <- c("EOMES", "FOXF1", "MIML1", "BMP4", "BMP7")
my_pseudotime_gene <- c("KRT19", "EOMES", "FABP1", "FABP2", "GATA4", "GSC", "FOXA1", "FOXA2", "SOX7", "SOX17", "HNF1B", "AFP", "CTNNB1", "GATA6", "GDF1", "GDF3", "HNF4A", "MIXL1", "SALL4")
my_pseudotime_gene <- c("BMP4", "CHRD", "FGF8", "FOXJ3", "GBX1", "NES", "NOG", "OTX2", "TP63", "PAX2", "PAX6", "SOX1", "TUBB3", "NCAM1", "VIM")
my_pseudotime_gene <- c("CDH1", "CXCL5", "DNMT3B", "HESX1", "IDO1", "LCK", "NANOG", "POU5F1", "SOX2", "TRIM22")
my_pseudotime_gene <- c("PAX6")
chupapergenes <- c("POU5F1", "NANOG", "SOX2", "DNMT3B", "NODAL", "EOMES", "ID1", "CDX1", "T", "MSX2", "CER1", "GATA4", "DKK4", "MYCT1", "POU2AF1", "PRDM1")
meso.markers <- c("T", "FOXC1", "GSC", "SNAI2", "SNAI1", "TBX6", "TWIST2", "NCLN", "BMP2", "BMP4", "BMP7", "INHBA", "NODAL", "TGFB1", "TGFB2", "CDH2", "GDF3", "NR5A2")
endo.markers <- c("CLDN6", "KRT19", "FOXA1", "FOXA2", "SOX7", "AFP", "CTNNB1", "SALL4", "CABP7", "CLDN1", "CPLX2", "HHEX", "LEFTY1", "PRDM1")
ecto.markers <- c("FGF8", "FOXJ3", "NES", "NOG", "OTX2", "PAX6", "SOX1", "TUBB3", "NCAM1", "VIM", "CDH9", "COL2A1", "DRD4", "LMX1A", "MAP2", "MYO3B", "NOS2", "NR2F1", "NR2F2", "PAPLN", "PAX3", "PRKCA", "SDC2", "TRPM8", "ZBTB16")
my_pseudotime_gene<- ecto.markers
my_top_genes <- row.names(subset(fData(cds_filtered), gene_short_name %in% my_pseudotime_gene))
cds_subset_top_genes <- cds_filtered[my_top_genes,]
cds_subset_top_genes <- orderCells(cds_subset_top_genes, reverse = FALSE)
plot_genes_in_pseudotime(cds_subset_top_genes, color_by = "Sample")


# plot expression by tSNE cluster for markers ----

cds <- clusterCells(cds, num_clusters = 5) #must change num_clusters per dataset. It has to be one less than the number of samples or clusters you expect. I don't know why.
plot_cell_clusters(cds, 1, 2, color = "Sample", markers = my_pseudotime_gene, cell_size = 0.5)
plot_cell_clusters(cds, 1, 2, color="Sample",cell_size = 0.5)

# BEAM------very slow-----
BEAM_res <- BEAM(cds, branch_point = 1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


# plot branched heatmap -----
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-50)),],
                            branch_point = 1,
                            num_clusters = 4, 
                            cores = 3,
                            use_gene_short_name = T,
                            show_rownames = T) # num_clusters here is arbitrarily set, different than num_clusters in clusterCells()

cds <- orderCells(cds, root_state = 1)

plot_multiple_branches_heatmap(cds[row.names(subset(BEAM_res,
                                                    qval < 1e-50)),],
                               branches = c(1, 2, 7, 6, 5),
                               branches_name = c("Days 0, 1, 2","Day 2", "Day 5", "Day 3", "Day 7"),
                               num_clusters = 2,
                               cores = 3,
                               use_gene_short_name = T,
                               show_rownames = T)
# cant take state "3" or "4"

# plot complex cell trajectory------
plot_complex_cell_trajectory(cds, color_by = 'Sample', show_branch_points = T,
                             cell_size = 1, cell_link_size = 1, root_states = c(1)) + scale_size(range = c(1, 1)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  theme (legend.position="left", legend.title=element_blank())

# plot cell trajectory with stat density-----
plot_cell_trajectory(cds, color_by = "Sample") + stat_density2d(color='black', h = 8, alpha=I(0.6), size=I(0.6))

plot_genes_in_pseudotime(cds, color_by = 'POU5F1')

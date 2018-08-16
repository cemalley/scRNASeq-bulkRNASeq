# Visualizing TSCAN pseudotime as an alternative to Monocle, and investigating if there are DE genes within sample age (day) by state.

library(TSCAN)
library(data.table)
library(Seurat)
library(ggplot2)

procdata <- preprocess(as.matrix(aggr6@raw.data))
lpsmclust <- exprmclust(as.matrix(aggr6@raw.data))
lpsorder <- TSCANorder(lpsmclust)
lpsorder
lpsmclust
TSCAN::plotmclust(lpsmclust, show_tree = T, show_cell_names = F)

orderTSCAN <- TSCAN::TSCANorder(lpsmclust, orderonly = F)
pseudotime_order_tscan <- as.character(orderTSCAN$sample_name)

orderTSCAN <- as.data.table(orderTSCAN)
orderTSCAN[,Day:=tstrsplit(sample_name, "\\.")[2]]
orderTSCAN

ggplot(data=orderTSCAN) + geom_point(aes(x=Pseudotime, y=State, color=Day), size=5) + labs(title="TSCAN Pseudotime for ddSeq day0-7")
ggplot(data=orderTSCAN) + geom_point(aes(x=Pseudotime, y=State, color=Day), size=5) + labs(title="TSCAN Pseudotime for ddSeq day0-7")


orderTSCAN.summary <- orderTSCAN[, .(.N), by = .(Day, State), ]
orderTSCAN.summary <- orderTSCAN.summary[order(Day)]
orderTSCAN.summary
sum(orderTSCAN.summary$N)

# within day3, between states 4 and 5----

day3.state4 <- unlist(orderTSCAN[Day=="day3" & State==4,c(sample_name)], use.names=F)
day3.state5 <- unlist(orderTSCAN[Day=="day3" & State==5,c(sample_name)], use.names=F)

aggr6.day3.state4 <- SubsetData(aggr6, cells.use=c(day3.state4))
aggr6.day3.state5 <- SubsetData(aggr6, cells.use=c(day3.state5))

aggr6.day3.state4 <- SetAllIdent(aggr6.day3.state4, id = "sample")
aggr6.day3.state4@meta.data$sample <- "State4"
aggr6.day3.state5 <- SetAllIdent(aggr6.day3.state5, id="sample")
aggr6.day3.state5@meta.data$sample <- "State5"

aggr6.day3.subset <- MergeSeurat(aggr6.day3.state4, aggr6.day3.state5)



aggr6.day3.subset <- NormalizeData(object = aggr6.day3.subset, normalization.method = "LogNormalize", scale.factor = 10000)
par(mfrow = c(1, 1))

aggr6.day3.subset <- FindVariableGenes(object = aggr6.day3.subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
aggr6.day3.subset <- ScaleData(object = aggr6.day3.subset, vars.to.regress = c("nUMI")) # regresses out nUMI

aggr6.day3.subset <- RunPCA(aggr6.day3.subset)
aggr6.day3.subset@ident
aggr6.day3.subset <- SetAllIdent(aggr6.day3.subset, id="sample")
PCAPlot(aggr6.day3.subset)
aggr6.day3.subset <- RunTSNE(aggr6.day3.subset)
TSNEPlot(aggr6.day3.subset)
anydiff.day3s <- FindMarkers(aggr6.day3.subset, ident.1 = 'State4')
anydiff.day3s
VlnPlot(aggr6.day3.subset, c('POU5F1'))
VlnPlot(aggr6.day3.subset, c('PAX6'))

VlnPlot(aggr6, ident.include = c('day5', 'day3'), c("PAX6"))
VlnPlot(aggr6, c("PAX6"))


# between state 6 and 7, all days----
state6 <- unlist(orderTSCAN[State==6,c(sample_name)], use.names=F)
state7 <- unlist(orderTSCAN[State==7,c(sample_name)], use.names=F)

aggr6.state6 <- SubsetData(aggr6, cells.use=c(state6))
aggr6.state7 <- SubsetData(aggr6, cells.use=c(state7))

aggr6.state6 <- SetAllIdent(aggr6.state6, id = "sample")
aggr6.state6@meta.data$sample <- "State6"
aggr6.state7 <- SetAllIdent(aggr6.state7, id="sample")
aggr6.state7@meta.data$sample <- "State7"

aggr6.subset <- MergeSeurat(aggr6.state6, aggr6.state7)

aggr6.subset <- NormalizeData(object = aggr6.subset, normalization.method = "LogNormalize", scale.factor = 10000)
par(mfrow = c(1, 1))

aggr6.subset <- FindVariableGenes(object = aggr6.subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
aggr6.subset <- ScaleData(object = aggr6.subset, vars.to.regress = c("nUMI")) # regresses out nUMI

aggr6.subset <- RunPCA(aggr6.subset)
aggr6.subset@ident
aggr6.subset <- SetAllIdent(aggr6.subset, id="sample")
PCAPlot(aggr6.subset)
aggr6.subset <- RunTSNE(aggr6.subset)
TSNEPlot(aggr6.subset)
anydiff <- FindMarkers(aggr6.subset, ident.1 = 'State6')
anydiff
VlnPlot(aggr6.subset, c('POU5F1'))
VlnPlot(aggr6.subset, c('LDHA'))
VlnPlot(aggr6.subset, c('PAX6'))


states4_7 <- MergeSeurat(aggr6.day3.state4, aggr6.day3.state5)
states4_7 <- MergeSeurat(states4_7, aggr6.state6)
states4_7 <- MergeSeurat(states4_7, aggr6.state7)

states4_7 <-NormalizeData(object = states4_7, normalization.method = "LogNormalize", scale.factor = 10000)
states4_7 <- FindVariableGenes(object = states4_7, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
states4_7 <- ScaleData(object = states4_7, vars.to.regress = c("nUMI")) # regresses out nUMI

states4_7 <- RunPCA(states4_7)
states4_7@ident
states4_7 <- SetAllIdent(states4_7, id="sample")
PCAPlot(states4_7)
states4_7 <- RunTSNE(states4_7)
TSNEPlot(states4_7)

VlnPlot(states4_7, c("PAX6"))
VlnPlot(states4_7, c('LDHA'))

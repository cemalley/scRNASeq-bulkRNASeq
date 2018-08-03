library(Seurat)
library(dplyr)
library(data.table)

aggr6 <- fread("Aggregation6.allLA_day0_day1_day2_day3_day5_day7.csv")
day0 <- as.data.frame(subset(aggr6, select=c(1, grep("day0", names(aggr6)))))
day1 <- as.data.frame(subset(aggr6, select=c(1, grep("day1", names(aggr6)))))
day2 <- as.data.frame(subset(aggr6, select=c(1, grep("day2", names(aggr6)))))
day3 <- as.data.frame(subset(aggr6, select=c(1, grep("day3", names(aggr6)))))
day5 <- as.data.frame(subset(aggr6, select=c(1, grep("day5", names(aggr6)))))
day7 <- as.data.frame(subset(aggr6, select=c(1, grep("day7", names(aggr6)))))

# reformat for seurat. this was necessary to retain metadata for each sample in the later, merged object with all days.
row.names(day0)<- day0$GeneId
day0 <- day0[-1]
day0 <- CreateSeuratObject(raw.data = day0, project = "NPC.day0", min.cells = 5, min.genes = 200)
day0@meta.data$sample <- "day0"
day0 <- NormalizeData(day0)
day0 <- ScaleData(day0, display.progress = F)

row.names(day1)<- day1$GeneId
day1 <- day1[-1]
day1 <- CreateSeuratObject(raw.data = day1, project = "NPC.day1", min.cells = 5, min.genes = 200)
day1@meta.data$sample <- "day1"
day1 <- NormalizeData(day1)
day1 <- ScaleData(day1, display.progress = F)

row.names(day2)<- day2$GeneId
day2 <- day2[-1]
day2 <- CreateSeuratObject(raw.data = day2, project = "NPC.day2", min.cells = 5, min.genes = 200)
day2@meta.data$sample <- "day2"
day2 <- NormalizeData(day2)
day2 <- ScaleData(day2, display.progress = F)

row.names(day3)<- day3$GeneId
day3 <- day3[-1]
day3 <- CreateSeuratObject(raw.data = day3, project = "NPC.day3", min.cells = 5, min.genes = 200)
day3@meta.data$sample <- "day3"
day3 <- NormalizeData(day3)
day3 <- ScaleData(day3, display.progress = F)

row.names(day5)<- day5$GeneId
day5 <- day5[-1]
day5 <- CreateSeuratObject(raw.data = day5, project = "NPC.day5", min.cells = 5, min.genes = 200)
day5@meta.data$sample <- "day5"
day5 <- NormalizeData(day5)
day5 <- ScaleData(day5, display.progress = F)

row.names(day7)<- day7$GeneId
day7 <- day7[-1]
day7 <- CreateSeuratObject(raw.data = day7, project = "NPC.day7", min.cells = 5, min.genes = 200)
day7@meta.data$sample <- "day7"
day7 <- NormalizeData(day7)
day7 <- ScaleData(day7, display.progress = F)

##

aggr6.seurat <- MergeSeurat(day0, day1, project="NPC", min.cells = 5, min.genes = 200)
aggr6.seurat <- MergeSeurat(aggr6.seurat, day2, project="NPC", min.cells = 5, min.genes = 200)
aggr6.seurat <- MergeSeurat(aggr6.seurat, day3, project="NPC", min.cells = 5, min.genes = 200)
aggr6.seurat <- MergeSeurat(aggr6.seurat, day5, project="NPC", min.cells = 5, min.genes = 200)
aggr6.seurat <- MergeSeurat(aggr6.seurat, day7, project="NPC", min.cells = 5, min.genes = 200)

VlnPlot(aggr6.seurat, c("nGene", "nUMI"), nCol = 2)

aggr6 <- aggr6.seurat

aggr6 <- SetAllIdent(object = aggr6, id = "sample")

#normalize
VlnPlot(object = aggr6, features.plot = c("nGene", "nUMI"), nCol = 2)

par(mfrow = c(1, 1))
GenePlot(object = aggr6, gene1 = "nUMI", gene2 = "nGene")

#mito genes - if you want to find and remove mitochondrial genes past a certain threshold of expression
mito.genes <- grep(pattern = "^MT-", x = rownames(x = aggr6@data), value = TRUE)
percent.mito <- Matrix::colSums(aggr6@raw.data[mito.genes, ])/Matrix::colSums(aggr6@raw.data)
aggr6 <- AddMetaData(object = aggr6, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = aggr6, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

aggr6 <- FilterCells(object = aggr6, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(7000))


aggr6 <- NormalizeData(object = aggr6, normalization.method = "LogNormalize", scale.factor = 10000)
par(mfrow = c(1, 1))

aggr6 <- FindVariableGenes(object = aggr6, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = aggr6@var.genes) #3682 genes for normalized aggregation 6

aggr6 <- ScaleData(object = aggr6, vars.to.regress = c("nUMI")) # regresses out nUMI

aggr6 <- RunPCA(object = aggr6, pc.genes = aggr6@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

VizPCA(object = aggr6, pcs.use = 1:2)

PCAPlot(object = aggr6, dim.1 = 1, dim.2 = 2)
PCAPlot(object = aggr6, dim.1 = 3, dim.2 = 4)
PCAPlot(object = aggr6, dim.1 = 5, dim.2 = 6)

## saved rdata: Figures/Seurat/Aggregation6.Normalized.RData

par(mfrow = c(2, 2))
PCHeatmap(object = aggr6, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = aggr6, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

aggr6 <- JackStraw(object = aggr6, num.replicate = 100)

JackStrawPlot(object = aggr6, PCs = 1:7)

PCElbowPlot(object = aggr6)

aggr6 <- FindClusters(object = aggr6, reduction.type = "pca", dims.use = 1:12, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = aggr6)

aggr6 <- SetAllIdent(object = aggr6, id = "sample")
aggr6 <- RunTSNE(object = aggr6, dims.use = 1:8, do.fast = TRUE)
tsneplotdata <- TSNEPlot(object = aggr6)

prop.table(x = table(aggr6@ident))

#

# cell cycle analysis----
g1.s.genes <- g1.s <- c("ACD", "ACYP1", "ADAMTS1", "ANKRD10", "APEX2", "ARGLU1", "ATAD2", "BARD1", "BRD7", "C1orf63", "C7orf41", "C14orf142", "CAPN7", "CASP2", "CASP8AP2", "CCNE1", "CCNE2", "CDC6", "CDC25A", "CDCA7", "CDCA7L", "CEP57", "CHAF1A", "CHAF1B", "CLSPN", "CREBZF", "CTSD", "DIS3", "DNAJC3", "DONSON", "DSCC1", "DTL", "E2F1", "EIF2A", "ESD", "FAM105B", "FAM122A", "FLAD1", "GINS2", "GINS3", "GMNN", "HELLS", "HOXB4", "HRAS", "HSF2", "INSR", "INTS8", "IVNS1ABP", "KIAA1147", "KIAA1586", "LNPEP", "LUC7L3", "MCM2", "MCM4", "MCM5", "MCM6", "MDM1", "MED31", "MRI1", "MSH2", "NASP", "NEAT1", "NKTR", "NPAT", "NUP43", "ORC1", "OSBPL6", "PANK2", "PCDH7", "PCNA", "PLCXD1", "PMS1", "PNN", "POLD3", "RAB23", "RECQL4", "RMI2", "RNF113A", "RNPC3", "SEC62", "SKP2", "SLBP", "SLC25A36", "SNHG10", "SRSF7", "SSR3", "TAF15", "TIPIN", "TOPBP1", "TRA2A", "TTC14", "UBR7", "UHRF1", "UNG", "USP53", "VPS72", "WDR76", "ZMYND19", "ZNF367", "ZRANB2")
s.genes <- c("ABCC5", "ABHD10", "ANKRD18A", "ASF1B", "ATAD2", "BBS2", "BIVM", "BLM", "BMI1", "BRCA1", "BRIP1", "C5orf42", "C11orf82", "CALD1", "CALM2", "CASP2", "CCDC14", "CCDC84", "CCDC150", "CDC7", "CDC45", "CDCA5", "CDKN2AIP", "CENPM", "CENPQ", "CERS6", "CHML", "COQ9", "CPNE8", "CREBZF", "CRLS1", "DCAF16", "DEPDC7", "DHFR", "DNA2", "DNAJB4", "DONSON", "DSCC1", "DYNC1LI2", "E2F8", "EIF4EBP2", "ENOSF1", "ESCO2", "EXO1", "EZH2", "FAM178A", "FANCA", "FANCI", "FEN1", "GCLM", "GOLGA8A", "GOLGA8B", "H1F0", "HELLS", "HIST1H2AC", "HIST1H4C", "INTS7", "KAT2A", "KAT2B", "KDELC1", "KIAA1598", "LMO4", "LYRM7", "MAN1A2", "MAP3K2", "MASTL", "MBD4", "MCM8", "MLF1IP", "MYCBP2", "NAB1", "NEAT1", "NFE2L2", "NRD1", "NSUN3", "NT5DC1", "NUP160", "OGT", "ORC3", "OSGIN2", "PHIP", "PHTF1", "PHTF2", "PKMYT1", "POLA1", "PRIM1", "PTAR1", "RAD18", "RAD51", "RAD51AP1", "RBBP8", "REEP1", "RFC2", "RHOBTB3", "RMI1", "RPA2", "RRM1", "RRM2", "RSRC2", "SAP30BP", "SLC38A2", "SP1", "SRSF5", "SVIP", "TOP2A", "TTC31", "TTLL7", "TYMS", "UBE2T", "UBL3", "USP1", "ZBED5", "ZWINT")
g2.m.genes<- c("ANLN", "AP3D1", "ARHGAP19", "ARL4A", "ARMC1", "ASXL1", "ATL2", "AURKB", "BCLAF1", "BORA", "BRD8", "BUB3", "C2orf69", "C14orf80", "CASP3", "CBX5", "CCDC107", "CCNA2", "CCNF", "CDC16", "CDC25C", "CDCA2", "CDCA3", "CDCA8", "CDK1", "CDKN1B", "CDKN2C", "CDR2", "CENPL", "CEP350", "CFD", "CFLAR", "CHEK2", "CKAP2", "CKAP2L", "CYTH2", "DCAF7", "DHX8", "DNAJB1", "ENTPD5", "ESPL1", "FADD", "FAM83D", "FAN1", "FANCD2", "G2E3", "GABPB1", "GAS1", "GAS2L3", "H2AFX", "HAUS8", "HINT3", "HIPK2", "HJURP", "HMGB2", "HN1", "HP1BP3", "HRSP12", "IFNAR1", "IQGAP3", "KATNA1", "KCTD9", "KDM4A", "KIAA1524", "KIF5B", "KIF11", "KIF20B", "KIF22", "KIF23", "KIFC1", "KLF6", "KPNA2", "LBR", "LIX1L", "LMNB1", "MAD2L1", "MALAT1", "MELK", "MGAT2", "MID1", "MIS18BP1", "MND1", "NCAPD3", "NCAPH", "NCOA5", "NDC80", "NEIL3", "NFIC", "NIPBL", "NMB", "NR3C1", "NUCKS1", "NUMA1", "NUSAP1", "PIF1", "PKNOX1", "POLQ", "PPP1R2", "PSMD11", "PSRC1", "RANGAP1", "RCCD1", "RDH11", "RNF141", "SAP30", "SKA3", "SMC4", "STAT1", "STIL", "STK17B", "SUCLG2", "TFAP2A", "TIMP1", "TMEM99", "TMPO", "TNPO2", "TOP2A", "TRAIP", "TRIM59", "TRMT2A", "TTF2", "TUBA1A", "TUBB", "TUBB2A", "TUBB4B", "TUBD1", "UACA", "UBE2C", "VPS25", "VTA1", "WSB1", "ZNF587", "ZNHIT2")
m.genes <- c("AHI1", "AKIRIN2", "ANKRD40", "ANLN", "ANP32B", "ANP32E", "ARHGAP19", "ARL6IP1", "ASXL1", "ATF7IP", "AURKA", "BIRC2", "BIRC5", "BUB1", "CADM1", "CCDC88A", "CCDC90B", "CCNA2", "CCNB2", "CDC20", "CDC25B", "CDC27", "CDC42EP1", "CDCA3", "CENPA", "CENPE", "CENPF", "CEP55", "CFLAR", "CIT", "CKAP2", "CKAP5", "CKS1B", "CKS2", "CNOT10", "CNTROB", "CTCF", "CTNNA1", "CTNND1", "DEPDC1", "DEPDC1B", "DIAPH3", "DLGAP5", "DNAJA1", "DNAJB1", "DR1", "DZIP3", "E2F5", "ECT2", "FAM64A", "FOXM1", "FYN", "G2E3", "GADD45A", "GAS2L3", "GOT1", "GRK6", "GTSE1", "HCFC1", "HMG20B", "HMGB3", "HMMR", "HN1", "HP1BP3", "HPS4", "HS2ST1", "HSPA8", "HSPA13", "INADL", "KIF2C", "KIF5B", "KIF14", "KIF20B", "KLF9", "LBR", "LMNA", "MCM4", "MDC1", "MIS18BP1", "MKI67", "MLLT4", "MZT1", "NCAPD2", "NCOA5", "NEK2", "NUF2", "NUP35", "NUP98", "NUSAP1", "ODF2", "ORAOV1", "PBK", "PCF11", "PLK1", "POC1A", "POM121", "PPP1R10", "PRPSAP1", "PRR11", "PSMG3", "PTP4A1", "PTPN9", "PWP1", "QRICH1", "RAD51C", "RANGAP1", "RBM8A", "RCAN1", "RERE", "RNF126", "RNF141", "RNPS1", "RRP1", "SEPHS1", "SETD8", "SFPQ", "SGOL2", "SHCBP1", "SMARCB1", "SMARCD1", "SPAG5", "SPTBN1", "SRF", "SRSF3", "SS18", "SUV420H1", "TACC3", "THRAP3", "TLE3", "TMEM138", "TNPO1", "TOMM34", "TPX2", "TRIP13", "TSG101", "TSN", "TTK", "TUBB4B", "TXNDC9", "TXNRD1", "UBE2D3", "USP13", "USP16", "VANGL1", "WIBG", "WSB1", "YWHAH", "ZC3HC1", "ZFX", "ZMYM1", "ZNF207")
m.g1.genes <- c("AGFG1", "AGPAT3", "AKAP13", "AMD1", "ANP32E", "ANTXR1", "BAG3", "BTBD3", "CBX3", "CDC42", "CDK7", "CDKN3", "CEP70", "CNIH4", "CTR9", "CWC15", "DCP1A", "DCTN6", "DEXI", "DKC1", "DNAJB6", "DSP", "DYNLL1", "EIF4E", "ELP3", "FAM60A", "FAM189B", "FOPNL", "FOXK2", "FXR1", "G3BP1", "GATA2", "GNB1", "GRPEL1", "GSPT1", "GTF3C4", "HIF1A", "HMG20B", "HMGCR", "HSD17B11", "HSPA8", "ILF2", "JMJD1C", "KDM5B", "KIAA0586", "KIF5B", "KPNB1", "KRAS", "LARP1", "LARP7", "LRIF1", "LYAR", "MORF4L2", "MRPL19", "MRPS2", "MRPS18B", "MSL1", "MTPN", "NCOA3", "NFIA", "NFIC", "NUCKS1", "NUFIP2", "NUP37", "ODF2", "OPN3", "PAK1IP1", "PBK", "PCF11", "PLIN3", "PPP2CA", "PPP2R2A", "PPP6R3", "PRC1", "PSEN1", "PTMS", "PTTG1", "RAD21", "RAN", "RHEB", "RPL13A", "SLC39A10", "SNUPN", "SRSF3", "STAG1", "SYNCRIP", "TAF9", "TCERG1", "TLE3", "TMEM138", "TOB2", "TOP1", "TROAP", "TSC22D1", "TULP4", "UBE2D3", "VANGL1", "VCL", "WIPF2", "WWC1", "YY1", "ZBTB7A", "ZCCHC10", "ZNF24", "ZNF281", "ZNF593")

aggr6 <- CellCycleScoring(object = aggr6, s.genes = s.genes, g2m.genes = g2.m.genes, set.ident = TRUE)
head(x =aggr6@meta.data)

ectoderm.genes <- c("VIM", "MAP2", "OTX2", "PAX6", "SOX1", "TUBB3", "COL2A1")
endoderm.genes <- c("CLDN6", "KRT19", "EOMES", "GSC", "FOXA1", "FOXA2", "SOX7", "AFP", "CTNNB1", "GDF3", "SALL4", "CABP7", "CLDN1", "CPLX2", "FOXA1", "FOXA2", "HHEX", "LEFTY1", "NODAL", "PRDM1")


#day0 de genes----
RidgePlot(object = aggr6, features.plot = c("ESRG", "FOXD3-AS1", "POU5F1", "LINC00678"), nCol = 2, group.by = "Phase")
RidgePlot(object = aggr6, features.plot = c(ectoderm.genes), nCol = 4, group.by = "sample")
RidgePlot(object = aggr6, features.plot = c(endoderm.genes), nCol = 8, group.by = "sample")
mesoderm.genes <- c("T", "EOMES", "FOXC1", "GSC", "SNAI2", "SNAI1", "TBX6", "TWIST2", "NCLN", "BMP2", "BMP4", "BMP6", "BMP7", "INHBA", "INHBB", "NODAL", "TGFB1", "TGFB2", "CDH2", "GDF3", "NR5A2")
RidgePlot(object = aggr6, features.plot = c(mesoderm.genes), nCol = 8, group.by = "sample")


#mesoderm de genes----
RidgePlot(object = aggr6, features.plot = c("DKK1", "TNFRSF11B", "MSX1", "NKD1"), nCol = 2, group.by = "Phase")
RidgePlot(object = aggr6, features.plot = c("DKK1", "TNFRSF11B", "MSX1", "NKD1"), nCol = 2, group.by = "sample")


# cell cycle scoring---
cell.cycle.scoring.out <- as.data.frame(aggr6@meta.data)
cell.cycle.scoring.out$barcode <- row.names(cell.cycle.scoring.out)
cell.cycle.scoring.out <- as.data.table(cell.cycle.scoring.out)
fwrite(cell.cycle.scoring.out, "R_scripts/Cell.cycle.scoring.Seurat.Aggr11.csv", sep="\t", col.names = T, row.names = F, quote=F)


# feature plots------
meso.markers <- c("T", "FOXC1", "GSC", "SNAI2", "SNAI1", "TBX6","TWIST2", "NCLN", "BMP2", "BMP4", "BMP7", "INHBA", "NODAL", "TGFB1", "TGFB2", "CDH2", "GDF3","NR5A2")
FeaturePlot(aggr6, c(meso.markers),cols.use = c("grey","blue"))

endo.markers <- c("CLDN6", "KRT19", "GSC", "FOXA1", "FOXA2", "SOX7", "AFP", "CTNNB1", "GDF3", "SALL4", "CABP7", "CLDN1", "CPLX2", "FOXA1", "FOXA2", "HHEX", "LEFTY1", "NODAL", "PRDM1")
FeaturePlot(aggr6, c(endo.markers),cols.use = c("grey","blue"))

ecto.markers <- c("BMP4", "FGF8", "FOXJ3", "NES", "NOG", "OTX2", "PAX6", "SOX1", "TUBB3", "NCAM1", "VIM", "CDH9", "COL2A1", "DRD4", "LMX1A", "MAP2", "MYO3B", "NOS2", "NR2F1", "NR2F2", "PAPLN", "PAX3", "PRKCA", "SDC2", "TRPM8", "ZBTB16")
FeaturePlot(aggr6, c(ecto.markers),cols.use = c("grey","blue"))

pluri.markers <- c("CDH1", "CXCL5", "DNMT3B", "HESX1", "IDO1", "LCK", "NANOG", "POU5F1", "SOX2", "TRIM22")

marker.heatmap <- DoHeatmap(aggr6, genes.use = c(meso.markers, endo.markers, ecto.markers), slim.col.label = TRUE, remove.key = TRUE, col.low = "green", col.high="red")
#marker.heatmap <- DoHeatmap(aggr6, genes.use = c(merged.genes), slim.col.label = TRUE, remove.key = TRUE, col.low = "green", col.high="red")

dim(FetchData(aggr6,c("ident",c(meso.markers, endo.markers, ecto.markers))))
length(c(meso.markers, endo.markers, ecto.markers))
heatmap.data.out <- FetchData(aggr6,c("ident",c(meso.markers, endo.markers, ecto.markers)))
heatmap.data.out$Barcode <- row.names(heatmap.data.out)
heatmap.data.out<- as.data.table(heatmap.data.out)
heatmap.data.out <- heatmap.data.out[,c(65, 1:64)]
fwrite(heatmap.data.out, "R_scripts/Seurat/Heatmap.data.csv", sep="\t", col.names = T, row.names = F, quote=F)


#grab genes used in whole seurat object
genelist <- row.names(as.data.frame(as.matrix(aggr6@data)))
length(grep("^MT-", genelist))
length(grep("^RP", genelist))
genelist <- as.data.table(unlist(genelist, use.names=F))
fwrite(genelist, "R_scripts/Seurat/Genelist.csv", sep="\t", col.names = F, row.names = F, quote=F)
# saved aggr6 seurat object R_scripts/Seurat/Aggr6.Seurat.Object.RData
#
lineages <- lineages[,c("Barcode", "Max")]

ordered.barcodes <- as.data.frame(unlist(c(row.names(aggr6@meta.data)), use.names=F))
names(ordered.barcodes) <- "Barcode"
ordered.barcodes$Barcode[1]
lineages$Barcode.ord <- factor(lineages$Barcode, levels=row.names(aggr6@meta.data))



# may 9 marker list-----
markers <- c("DPYSL5", "PRTG", "GREB1L", "TPBG", "CRB2", "SKI", "SPON1", "SIX3", "WLS", "GPC3")

# may 10 marker list ----
markers <-c("CPT1A", "IRS1", "PDK1", "INSR", "PRKAG1", "PPP2R5A", "PIK3R2", "PIK3R1", "HMGCR",
            "FOXO1", "IGF1R", "PPP2R2B", "AKT2", "MLYCD", "CPT1A", "ACSL1", "ACSL4", "ACOX1", "ACAA1",
            "ACAA1", "ACAA2", "ACAT1", "ACAT2", "ACAD9", "ACAD10", "ACADM",
            "ACADS", "ACADSB", "ACADVL", "EHHADH", "GCDH", "ACOX1", "ACOX2", "ACOX3", "ACSBG1",
            "ACSL1", "ACSL3", "ACSL4", "ACSL6", "ACSM3", "ACSM4",
            "ACOT7", "ACOT8", "ACOT9", "CPT1A", "CPT1C", "CPT2", "CRAT", "CROT", "ALDH2",
            "DECR1", "DECR2", "ECHS1", "ECI2", "HADHA", "MCEE", "MUT", "PECR", "PPA1",
            "FABP3", "FABP5", "FABP6", "FABP7", "SLC27A1", "SLC27A2", "SLC27A3", "SLC27A4",
            "SLC27A5", "SLC27A6", "FASN", "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKACA", "PRKACB",
            "PRKAG1", "PRKAG2", "BDH1", "BDH2", "HMGCL", "HMGCS1", "GK",
            "LIPE", "LPL", "GAPDH", "ACTB", "B2M", "RPLP0", "RPL13A", "TUBB", "CFL1")
# heatmap----
DoHeatmap(aggr6, genes.use=c(markers), col.low = "green", col.high="red", remove.key = T, slim.col.label = TRUE)

#violin plots----
VlnPlot(aggr6, c(markers))


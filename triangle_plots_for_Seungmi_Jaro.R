library(ggtern)
library(ggplot2)
library(data.table)
library(readxl)
library(Seurat)
library(ggthemes)

# Seungmi data ------

# using old list
marker.lists <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/Trilineage_markers.xlsx'))
meso.genes <- unique(na.omit(marker.lists$meso))
endo.genes <- unique(na.omit(marker.lists$endo))
ecto.genes <- unique(na.omit(marker.lists$ecto))

load("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData")

names(is019@meta.data)
unique(is019@meta.data$clusters)

Idents(is019) <- 'clusters'

cluster.averages <- AverageExpression(is019)

data <- cluster.averages$RNA
data$gene <- row.names(data)
data <- subset(data, select=c(8, 1:7))
data <- as.data.table(data)
data

meso.dt <- data[gene %in% meso.genes,]
endo.dt <- data[gene %in% endo.genes,]
ecto.dt <- data[gene %in% ecto.genes,]

meso.means <- colMeans(meso.dt[,-1])
endo.means <- colMeans(endo.dt[,-1])
ecto.means <- colMeans(ecto.dt[,-1])

meso.means


data.tern <- rbind(meso.means, endo.means, ecto.means)
data.tern <- as.data.table(t(data.tern))
data.tern
names(data.tern) <- c('Mesoderm','Endoderm','Ectoderm')



data.tern[,Cluster:= c(names(data)[2:8])]

data.tern.avgs <- data.tern

ggtern(data.tern, aes(Mesoderm, Endoderm, Ectoderm, color=Cluster)) + 
  theme_hc() +
  geom_point(size=3)+theme_rgbw() +
  theme_grid() +  
  geom_Lmark(colour='darkblue') + 
  geom_Tmark(colour='darkred') +
  geom_Rmark(colour='darkgreen') 

# plot all cells-----

data <- as.data.frame(as.matrix(is019@assays$RNA@data))
data[1:10,1:10]
data$gene <- row.names(data)
data <- subset(data, select=c(8951, 1:8950))
data <- as.data.table(data)

meso.dt <- data[gene %in% meso.genes,]
endo.dt <- data[gene %in% endo.genes,]
ecto.dt <- data[gene %in% ecto.genes,]
meso.means <- colMeans(meso.dt[,-1])
endo.means <- colMeans(endo.dt[,-1])
ecto.means <- colMeans(ecto.dt[,-1])

data.tern <- rbind(meso.means, endo.means, ecto.means)
data.tern <- as.data.table(t(data.tern))
data.tern
names(data.tern) <- c('Mesoderm','Endoderm','Ectoderm')


data.tern[,Cluster:= mergedSeuratProcessed@meta.data$orig.ident]

data.tern.cells <- data.tern
data.tern.cells[,Type:='cell']

# averaging.
cluster.averages <- AverageExpression(mergedSeuratProcessed)

data <- cluster.averages$RNA
data$gene <- row.names(data)
data <- subset(data, select=c(7, 1:6))
data <- as.data.table(data)
data

meso.dt <- data[gene %in% meso.genes,]
endo.dt <- data[gene %in% endo.genes,]
ecto.dt <- data[gene %in% ecto.genes,]

meso.means <- colMeans(meso.dt[,-1])
endo.means <- colMeans(endo.dt[,-1])
ecto.means <- colMeans(ecto.dt[,-1])

meso.means


data.tern <- rbind(meso.means, endo.means, ecto.means)
data.tern <- as.data.table(t(data.tern))
data.tern
names(data.tern) <- c('Mesoderm','Endoderm','Ectoderm')



data.tern[,Sample:= c(names(data)[2:7])]

data.tern.avgs <- data.tern
data.tern.avgs[,Type := 'avg']

data.tern.cells <- rbind(data.tern.cells, data.tern.avgs)
data.tern.cells[,size:=c(rep(1, 12744), rep(2, 6))]
data.tern.cells[,alpha:=c(rep(2, 12744), rep(2.1, 6))]
data.tern.cells

ggtern(data.tern.cells, aes(Mesoderm, Endoderm, Ectoderm, color=Sample, size=size)) + 
  theme_bw() +
  geom_point()+
  guides(size=FALSE)

ggtern(data.tern.cells[Type=='avg',], aes(Mesoderm, Endoderm, Ectoderm, color=Sample, size=size)) + 
  theme_bw() +
  geom_point()+
  guides(size=FALSE)

# added black stroke on avg points in illustrator.


# Jaro--------------------------
load("/Volumes/ncatssctl/NGS_related/Chromium/IS021/Ben/seurat_processed_cell-labels_10082019.RData")

trophectoderm <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/placental genes.txt', header=F)
trophectoderm

Idents(mergedSeuratProcessed) <- 'orig.ident'
mergedSeuratProcessed@active.ident


data <- as.data.frame(as.matrix(mergedSeuratProcessed@assays$RNA@data))
data[1:10,1:10]
data$gene <- row.names(data)
data <- subset(data, select=c(12745, 1:12744))
data <- as.data.table(data)

meso.dt <- data[gene %in% meso.genes,]
endo.dt <- data[gene %in% endo.genes,]
ecto.dt <- data[gene %in% ecto.genes,]
meso.means <- colMeans(meso.dt[,-1])
endo.means <- colMeans(endo.dt[,-1])
ecto.means <- colMeans(ecto.dt[,-1])

data.tern <- rbind(meso.means, endo.means, ecto.means)
data.tern <- as.data.table(t(data.tern))
data.tern
names(data.tern) <- c('Mesoderm','Endoderm','Ectoderm')


data.tern[,Sample:= mergedSeuratProcessed@meta.data$orig.ident]

data.tern.cells <- data.tern
data.tern.cells[,Type:='cell']

# re-ran averaging.
data.tern.avgs[,Type := 'avg']

data.tern.cells <- rbind(data.tern.cells, data.tern.avgs)
data.tern.cells[,size:=c(rep(1, 12744), rep(2, 6))]
data.tern.cells[,alpha:=c(rep(2, 12744), rep(2.1, 6))]
data.tern.cells

ggtern(data.tern.cells, aes(Mesoderm, Endoderm, Ectoderm, color=Sample, size=size)) + 
  theme_bw() +
  geom_point()+
  guides(size=FALSE)





#
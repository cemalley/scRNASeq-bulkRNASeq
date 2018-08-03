library(FateID) # see https://github.com/dgrun/FateID
library(Seurat)
library(data.table)

aggr6 <- fread("Aggregation6.allLA_day0_day1_day2_day3_day5_day7.csv")
aggr6 <- aggr6[!(grep('^RP', GeneId)),] #remove ribosomal proteins

x <- as.data.frame(aggr6)
row.names(x)<- aggr6$GeneId
x <- x[,-1]

# start FateID with new clusters defined by strong gene markers-----
pluri.markers <- c("ESRG", "NANOG", "POU5F1", "DUSP6", "ZFP42", "PODXL")
ecto.markers <- c("COL2A1", "PAX6", "SOX1", "NR2F2", "VIM")
meso.markers <- c("T", "HHEX", "NCLN")
endo.markers <- c("CTNNB1")
endpoint.markers <- c("PAX6", "COL2A1", "OTX2")

FMarker <- list(c(pluri.markers), c(meso.markers), c(endo.markers), c(ecto.markers))
xf <- getPart(x,FMarker,fthr=NULL,n=20)
head(xf$part)
head(xf$tar)
tar <- xf$tar
y <- xf$part
rc <- reclassify(x, y, tar, clthr=.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL, q=0.9)
y  <- rc$part
x  <- rc$xf

x <- getFeat(x,y,tar,fpv=0.01)

fb  <- fateBias(x, y, tar, z=NULL, minnr=10, minnrh=20, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)

dr  <- compdr(x, z=NULL, m=c("tsne","cmd","dm","lle"), k=c(2,3), lle.n=30, dm.sigma="local",
              dm.distance="euclidean", tsne.perplexity=30, seed=12345)

plotFateMap(y,dr,k=2,m="tsne")

plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t2", n="Pluripotency markers bias") #pluri markers
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t3", n="Mesoderm markers bias") # meso markers
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t4", n="Endoderm markers bias") # endo markers
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="t5", n="Ectoderm markers bias") # ecto markers
plotFateMap(y,dr,k=3,m="tsne",fb=fb,g="t5", n="Ectoderm markers bias") #3D
plotFateMap(y,dr,k=2,m="tsne",fb=fb,g="E", n="Fate bias entropy")

dev.off()
pr <- plotFateMap(y,dr,k=2,m="tsne",trthr=.33,fb=fb,prc=TRUE)
k <- impGenes(fb,"t5")
impGenes(fb,"t5")
impGenes(fb,"t4")
impGenes(fb,"t3")
impgenes2<- impGenes(fb,"t2")
row.names(impgenes2$d)
DoHeatmap(aggr6, genes.use=c(row.names(impgenes2$d)), col.low=c("green"), col.high=c("red"), slim.col.label = TRUE)


# partition by day
part.output <- as.data.frame(xf$part)
part.output$cell <- row.names(part.output)
part.output <- as.data.table(part.output)
part.output[,day:= tstrsplit(cell, "\\.")[2]]
part.output
library(ggplot2)
library(ggthemes)
ggplot(part.output) + geom_jitter(aes(x=part.output$day, y=part.output$`xf$part`), height = 0, width=0.1) + scale_y_continuous(breaks=c(1:5)) + theme_bw()+ labs(x="Day", y="Partition", title="RF partitions by day")


# find fate bias towards a particular lineage -- self-organizing map (SOM) of the pseudo-temporal expression profiles is computed
fs  <- filterset(x,n=pr$trc[["t5"]],minexpr=2,minnumber=1)
s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
set.seed(111111)
fcol <- sample(rainbow(max(y)))
n <-pr$trc[["t3"]]
plotheatmap(ps$nodes.z, xpart=y[n], xcol=fcol, ypart=unique(ps$nodes), xgrid=T, ygrid=TRUE, xlab=FALSE)
plotheatmap(ps$all.z, xpart=y[n], xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)
plotheatmap(ps$all.e, xpart=y[n], xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

g <- names(ps$nodes)[ps$nodes == 1]
plotexpression(fs, y, g, n, k=25, col=fcol, name="Node 1", cluster=FALSE, locreg=TRUE, alpha=.5, types=NULL)
dev.off()

# DE ------
thr <- .5
a   <- "t2"
b   <- "t3"
cl  <- c(2,3,4,5)
A <- rownames(fb$probs)[fb$probs[,a] > thr]
A <- A[y[A] %in% cl]
B <- rownames(fb$probs)[fb$probs[,b] > thr]
B <- B[y[B] %in% cl]

v <- log2(x+1)

de <- diffexpnb(v,A=A,B=B,DESeq=FALSE,norm=FALSE,vfit=NULL,locreg=FALSE)
de$res
dev.off()
plotdiffgenesnb(de,mthr=-4,lthr=0,Aname=a,Bname=b,padj=FALSE)
dev.off()
gene2gene(v,y,"PAX6","NANOG", plotnum=FALSE, fb=fb, tn="t3")
k <- impGenes(fb,"t2")

# plotting partitions----

melt.y <- as.data.frame(melt(y))
melt.y$cell <- row.names(melt.y)
melt.y <- as.data.table(melt.y)
melt.y[,day:=tstrsplit(cell, "\\.")[2]]

ggplot(data=melt.y) + geom_violin(aes(x=day, y=value))


# heatmap from fateid importance value genes---
load("Aggr6.Seurat.Object.RData")
genes <- c("ESRG", "CKB", "CDH1", "PODXL", "CD24", "HSPA8", "TDGF1", "EZR", "EIF5A", "HSPD1", "PSPIP1", "PRTG","MALAT1", "UCHL1", "KARS","EIF1","GLO1", "SRSF1","UGP2", "COX5A","MAD2L2","FOXD3-AS1","CD9","DPPA4",
           "CYCS","SLC16A1","GSTP1","SCNN1A","RARRES2","LECT1","TPBG","HSP90AB1","PTGES3","ACTG1","PPP1CC")

DoHeatmap(aggr6, genes.use=c(genes), col.low="#3660a5", col.mid = "white", col.high="#cd1f20", slim.col.label = TRUE)

# expression plots----
n <- pr$trc[["t5"]]
trc <- dptTraj(x,y,fb,trthr=.25,distance="euclidean",sigma=1000)
v
fs  <- filterset(v,n=n,minexpr=2,minnumber=1)
plotexpression(fs, y, g="PAX6", n=n, k=25, cluster=FALSE, locreg=TRUE, alpha=.5, types=substring(n, 20, 24))


# heatmap of z-scores---
set.seed(111111)
fcol <- sample(rainbow(max(y)))
s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
plotheatmap(ps$nodes.z, xpart=y[n], xcol=fcol, ypart=unique(ps$nodes), xgrid=FALSE, ygrid=TRUE, xlab=FALSE)

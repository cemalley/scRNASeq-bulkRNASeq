library(data.table)
library(ggthemes)
library(ggrepel)

# volcano plot, CMT2A21----
#res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB019/Analysis/DE/DE.H9noc_D56.vs.CMT2A21.csv')
res <- fread('~/Downloads/CMT2A_plus_D28noci/DE/All_DE_results/DE.H9noc_D28.vs.CMT2A21.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'CMT2A21'

# set p to 1*10^-300 where padj=0 for plotting gene labels more easily.
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

right <- subset(restoplot, ((log2FoldChange > 1) & (-log10(padj_reduce) > 5))) #0
left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj_reduce) > 5)))) #0

if(nrow(right)>15){
  right <- right[order(-log2FoldChange)]
  right <- right[1:15,]
}

if(nrow(left)>15){
  left <- left[order(log2FoldChange)]
  left <- left[1:15,]
}

restoplot[,threshold:=ifelse((-log10(padj_reduce) > 50) & log2FoldChange > 1, '#007F00', 
                             ifelse( (-log10(padj_reduce) > 50) & log2FoldChange < -1, 'red', 'gray' ))]

plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                            color=restoplot$threshold)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj_reduce), label=GeneId),
                  segment.size = 0.5)+
  geom_vline(aes(xintercept=1),color='brown')+
  geom_vline(aes(xintercept=-1),color='brown')+
  geom_hline(aes(yintercept=50), color='brown')+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: 266, DOWN: 258'), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  xlim(-8,8)+ylim(0,14) # may need to expand to -20,20 if the log2fold changes are that big
plot

# without lines or up/down-----
res <- fread('~/Downloads/CMT2A_plus_D28noci/DE/All_DE_results/DE.H9noc_D28.vs.CMT2A21.csv')
restoplot <- na.omit(res)
condition1 <- 'H9noc_D28'
condition2 <- 'CMT2A21'
res.sig <- fread('~/Downloads/CMT2A_plus_D28noci/DE/Significant_DE/DE.sig.H9noc_D28.vs.CMT2A21.csv')
res.sig
right <- res.sig[log2FoldChange>0,]
left <- res.sig[log2FoldChange<0,]
restoplot[,threshold:=ifelse(GeneId %in% left$GeneId, '#007F00', 
    ifelse( GeneId %in% right$GeneId, 'red', 'gray' ))]
restoplot[,padj_reduce:=padj]
restoplot[padj_reduce ==0,padj_reduce := (1*10^(-300))]

n.up <- nrow(right)
n.down <- nrow(left)

ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj_reduce)),
                                    color=restoplot$threshold)+
  geom_text(label="Higher in CMT2A21", x=6, y=0, color='red', vjust=1)+
  geom_text(label="Higher in H9noc D28", x=-6, y=0, color='#007F00', vjust=1)+
  geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj), label=GeneId),
                  segment.size = 0.5)+
  geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj), label=GeneId),
                  segment.size = 0.5)+
  labs(title=paste0('Volcano plot: ', condition1, ' vs ' ,
                    condition2, '. UP: ',n.up,'. DOWN: ',n.down), x='Effect size: log2(fold-change)',
       y='-log10(adjusted p-value)')+
  theme_bw()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  xlim(-8,8)+ylim(0,14)

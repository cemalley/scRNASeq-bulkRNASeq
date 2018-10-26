#syntax
#DESeq2_pipeline <- function(dds, condition1, condition2, testdir, volcanodir)

DESeq2_pipeline <- function(dds, condition1, condition2, testdir, volcanodir){
  library(data.table)
  library(DESeq2)
  library(ggrepel)
  library(ggplot2)
  
  res <- lfcShrink(dds, contrast=c("condition",condition2, condition1)) # default padjust method is Benjamini-Hochberg. must specify "later timepoint" versus "early timepoint" to be in the same log fold change sign direction as "early versus later" in results.
  
  reformat.res <- function(x, condition1, condition2){
    x <- as.data.frame(x)
    x$GeneId <- row.names(x)
    x <- as.data.table(x)
    x[,baseMean:=NULL]
    x[,baseMean1 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition1]))]
    x[,baseMean2 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition2]))]
    x <- na.omit(x)
    x <- x[order(rank(padj))]
    x <- x[,c("GeneId", "baseMean1", "baseMean2", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    names(x)[2:3] <- c(condition1, condition2)
    return(x)
  }
  
  res <- reformat.res(res, condition1, condition2)
  
  fwrite(x = res, file = paste0(testdir,'/','DE.', condition1, '.vs.', condition2,'.csv'))

  restoplot <- na.omit(res)
  right <- subset(restoplot, ((log2FoldChange > 1) & (-log10(padj) > 50)))
  left <- subset(restoplot, ( ((log2FoldChange < -1) & (-log10(padj) > 50))))

  plot <- ggplot(data=restoplot) + geom_point(aes(x=log2FoldChange, y=-log10(padj)),color='skyblue')+
    geom_text_repel(data=right, aes(x=log2FoldChange, y=-log10(padj), label=GeneId),
                    segment.size = 0.5,
                    nudge_x = 5 + right$log2FoldChange )+
    geom_text_repel(data=left, aes(x=log2FoldChange, y=-log10(padj), label=GeneId),
                    segment.size = 0.5,
                    nudge_x = -20 - left$log2FoldChange)+
    geom_vline(aes(xintercept=1),color='brown')+
    geom_vline(aes(xintercept=-1),color='brown')+
    geom_hline(aes(yintercept=50), color='brown')+
    labs(title=paste0('Volcano plot: ', condition1, ' vs ' , condition2), x='Effect size: log2(fold-change)', y='-log10(adjusted p-value)')+
    theme_bw()+
    xlim(-12,12)

  ggsave(plot, device='pdf', filename=paste0(condition1,'.vs.', condition2, 'volcano.pdf'), path=volcanodir, width=10, height=10, units='in')

}


# planning all DE tests----

celllines <- c('H9', 'Lonza')

times <- c('24hr', '72hr', 'Thaw')

conditions <- c('CE', 'CEPT', 'Y27', 'EDTA', 'CH')
expand.grid(celllines,conditions, times)


#i.e. H9_CE_24hr vs H9_CEPT_24hr


text <- capture.output( for (i in times){
  for (b in celllines){
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[2], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[3], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[1], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[3], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[2], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[3], '_', i, ' vs ', b, '_', conditions[4], '_', i, '\n'))
  cat(paste0(b, '_', conditions[3], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  cat(paste0(b, '_', conditions[4], '_', i, ' vs ', b, '_', conditions[5], '_', i, '\n'))
  }
})



tests_table <- data.table('test'=text)
tests_table[,c('condition1', 'condition2'):=tstrsplit(test, ' vs ')]

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Volcano'

cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]

  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE'
  volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Volcano'

  print(paste0('Now testing: ', condition1, ' vs ', condition2))
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)

  print(paste0('Finished: ', condition1, ' vs ', condition2))
}

stopCluster(cl)

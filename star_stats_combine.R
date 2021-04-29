# combine STAR alignment stats output files

library(data.table)
library(stringr)

#samples <- c("A04","B04","C04","D04")

#samples <- c('VJ3736_1_S16', 'VJ3736_2_S17', 'VJ3736_3_S18', 'VJ3736_4_S19', 'VJ3736_5_S20', 'VJ3736_6_S21')

#samples <- c('ISB008-A01_S1', 'ISB008-A02_S4', 'ISB008-C02_S6', 'ISB008-E03_S12')

#setwd('/Volumes/ncatssctl/Claire/data_requests/STAR_stats')
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis/STAR_stats/')
samples <- Sys.glob('*final*out')
samples <- gsub('_hg38Log.final.out','', samples)

stats <- data.table()

for (sample in samples){
  if (nrow(stats) > 0){
    stats.temp <- fread(paste0(sample, "_hg38Log.final.out"), header=F, sep="|", fill=T)
    stats.temp$V2 <- lapply(stats.temp$V2, function(x){gsub('\t','',x)})
    stats.temp <- as.data.table(t(stats.temp))
    names(stats.temp)<- c(unlist(stats.temp[1,]), use.names=F)
    stats.temp <- stats.temp[-1,]
    stats.temp <-stats.temp[,-c(5,8,23,28,32)]
    stats.temp$sample <- sample
    stats <- rbind(stats, stats.temp)
  }
  
  if(nrow(stats) < 1){
    stats <- fread(paste0(sample, "_hg38Log.final.out"), header=F, sep="|", fill=T)
    stats$V2 <- lapply(stats$V2, function(x){gsub('\t','',x)})
    stats <- as.data.table(t(stats))
    names(stats)<- c(unlist(stats[1,]), use.names=F)
    stats <- stats[-1,]
    stats <- stats[,-c(5,8,23,28,32)]
    stats$sample <- sample
  }
}
stats

fwrite(stats, 'ISB025_CT_STAR_stats.csv')


## combining across multiple runs and sequencing cores-----

#setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/quality_across_runs')
# 
# samples <- c('VJ3736_1_S16', 'VJ3736_2_S17', 'VJ3736_3_S18', 'VJ3736_4_S19', 'VJ3736_5_S20', 'VJ3736_6_S21', 'IS007_A04', 'IS007_B04', 'IS007_C04', 'IS007_D04',
#              'Tao_Jaro_PSC_16_D3_A1_1_S14', 'Tao_nociceptor_1_H9noc_D0_1', 'Robin_Sample_1', 'ISB003_30_S30')
# 
# stats <- data.table()
# 
# for (sample in samples){
#   if (nrow(stats) > 0){
#     stats.temp <- fread(paste0(sample, "_hg38Log.final.out"), header=F, sep="|", fill=T)
#     stats.temp$V2 <- lapply(stats.temp$V2, function(x){gsub('\t','',x)})
#     stats.temp <- as.data.table(t(stats.temp))
#     names(stats.temp)<- c(unlist(stats.temp[1,]), use.names=F)
#     stats.temp <- stats.temp[-1,]
#     stats.temp <-stats.temp[,-c(5,8,23,28,32)]
#     stats.temp$sample <- sample
#     stats <- rbind(stats, stats.temp)
#   }
#   
#   if(nrow(stats) < 1){
#     stats <- fread(paste0(sample, "_hg38Log.final.out"), header=F, sep="|", fill=T)
#     stats$V2 <- lapply(stats$V2, function(x){gsub('\t','',x)})
#     stats <- as.data.table(t(stats))
#     names(stats)<- c(unlist(stats[1,]), use.names=F)
#     stats <- stats[-1,]
#     stats <- stats[,-c(5,8,23,28,32)]
#     stats$sample <- sample
#   }
# }



stats$`Uniquely mapped reads %` <- lapply(stats$`Uniquely mapped reads %`, function(x){gsub('%','',x)})
stats


library(ggplot2)
attach(stats)
library(scales)
library(ggthemes)


stats[,mapped_reads := as.numeric(`Uniquely mapped reads %`)]



#ggplot(data=stats) + geom_bar(aes(y=`Total reads, average per sample`, x=`Run`, fill=`Platform`), stat = "identity") + scale_y_continuous(name="Total reads, average per sample", labels = comma)
#+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  geom_text(aes(x=`Run`, y=`Total reads, average per sample`, label = `Cell type`), vjust = -0.5)

# ggplot(data=stats, aes(x=sample, y=`Uniquely mapped reads %`)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   geom_text(aes(x=sample, y=`Uniquely mapped reads %`, label = formatC(`Uniquely mapped reads %`, digits=0, format="d") ), vjust = -0.5)+
#   labs(x='Sample', title='GTEx uniquely mapped reads %')

stats_label <- stats[mapped_reads < 70,]

library(ggrepel)

p <- ggplot() + geom_point(data=stats, aes(x=reorder(sample,-mapped_reads), y=mapped_reads)) +theme_clean()+ theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15), axis.text.y = element_text(size=12)) + 
  geom_text_repel(data=stats_label, aes(x=sample, y=mapped_reads, label = formatC(mapped_reads, digits=0, format="d"))) +
  #geom_text(data=stats_label, aes(x=sample, y=mapped_reads, label = formatC(mapped_reads, digits=0, format="d"), vjust=-0.5 ))+
  labs(x='Sample', title='GTEx mapping quality', y='Uniquely mapped reads %')+ylim(0,100)

p

library(grid)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

cat(stats[mapped_reads < 70,sample], sep='\n')

# facilities <- readxl::read_xlsx('facilities.xlsx')
# 
# stats <- merge(stats, facilities, by='sample')
# #stats <- as.data.table(stats)
# #stats[sample.ordered:= factor(sample, levels=c("IS007_A04", "IS007_B04", "IS007_C04", "IS007_D04", "Robin_Sample_1", "ISB003_30_S30", "Tao_Jaro_PSC_16_D3_A1_1_S14", "Tao_nociceptor_1_H9noc_D0_1", "VJ3736_1_S16", "VJ3736_2_S17", "VJ3736_3_S18", "VJ3736_4_S19", "VJ3736_5_S20", "VJ3736_6_S21"))]
# 
# attach(stats)
# library(forcats)
# 
# stats %>% mutate(sample = fct_relevel(sample, "IS007_A04", "IS007_B04", "IS007_C04", "IS007_D04", "Robin_Sample_1", "ISB003_30_S30", "Tao_Jaro_PSC_16_D3_A1_1_S14", "Tao_nociceptor_1_H9noc_D0_1", "VJ3736_1_S16", "VJ3736_2_S17", "VJ3736_3_S18", "VJ3736_4_S19", "VJ3736_5_S20", "VJ3736_6_S21")) %>%
#   ggplot(aes(x=sample, y=`Uniquely mapped reads %`)) + geom_bar(stat="identity", aes(fill=stats$Facility)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   geom_text(aes(x=sample, y=`Uniquely mapped reads %`, label = formatC(`Uniquely mapped reads %`, digits=0, format="d") ), vjust = -0.5)+
#   labs(x='Sample',title='Alignment quality across bulk RNASeq runs')+
#   scale_fill_discrete(name='Facility')
# 
# 
# fwrite(stats, 'STAR_alignment_quality_across_runs.csv')
#
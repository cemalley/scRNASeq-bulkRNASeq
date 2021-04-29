# merge all DE test results, nosiRNA day0 vs. all----
# I did AN2.DayX vs Sample.DayX.

setwd('/Volumes/ncatssctl/NGS_related/RASLseq/RASLResults/R_analysis/ISRL004-5_DE/nosiRNA/day0')
files <- Sys.glob('DE.nosiRNA.Control.Day0*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)

condition1 <- str_split_fixed(files[1], '\\.', Inf)[2]
condition2 <- str_split_fixed(files[1], '\\.', Inf)[7]
day <- str_split_fixed(files[1], '\\.', Inf)[9]
run <- str_split_fixed(files[1], '\\.', Inf)[5]
siRNA <- paste(str_split_fixed(files[1],'\\.',Inf)[7],str_split_fixed(files[1],'\\.',Inf)[8], sep='-')
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
data.merged[,day:= paste0(day)]
data.merged[,run:= paste0(run)]
data.merged[,siRNA:= paste0(siRNA)]
names(data.merged)[2] <- 'nosiRNA_mean'
names(data.merged)[3] <- 'Alt_mean'
data.merged

library(stringr)

for (file in files[2:118]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[7]
  day <- str_split_fixed(file, '\\.', Inf)[9]
  run <- str_split_fixed(file, '\\.', Inf)[5]
  siRNA <- paste(str_split_fixed(file[1],'\\.',Inf)[7],str_split_fixed(file[1],'\\.',Inf)[8], sep='-')
  data[,test:= paste0(condition1, '.vs.', condition2)]
  data[,day:= paste0(day)]
  data[,run:= paste0(run)]
  data[,siRNA:= paste0(siRNA)]
  names(data)[2] <- 'nosiRNA_mean'
  names(data)[3] <- 'Alt_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
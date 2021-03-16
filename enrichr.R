library(enrichR, quietly = T)
library(data.table)
# Written first by John Braisted, edited by Claire Malley

combineResults <-function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                        oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
      names(currdf) <- col_names
      fullRes <- rbind(fullRes, currdf)
    }
  }
  return(fullRes)
}
getSigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                           oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      if(res[j,4]<0.1) {				
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        fullSigRes <- rbind(fullSigRes, currdf)
      }
    }
  }
  return(fullSigRes)
}
getTop10SigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  topSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                          oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  for (i in 1:length(libs)) {
    #loop over result, append lib name to front, save only adjP <0.1
    res<-data.frame(enr[i])
    rowCnt<-dim(res)[2]
    # only consider the top 10 scores, only take adjP < 0.1
    for(j in 1:min(10,rowCnt)) {
      print(paste("**",res[j,4],"**"))
      if((!is.na(res[j,4])) & (res[j,4]<0.1)) {
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        topSigRes <- rbind(topSigRes, currdf)
      }
    }
  }
  return(topSigRes)
}
unpivotGeneLists <- function(fsr) {
  geneFrame<- data.frame(libName=character(), term=character(), adjP=double(), comb_score=double(), gene_symbol=character(), stringsAsFactors=FALSE)
  currRow <- 1
  mainRow <- 1
  for(r in 1:(dim(fsr)[1])){
    genes <- unlist(strsplit(as.character(fsr[r,"gene_list"]),';'))		
    for(gene in 1:length(genes)){
      geneFrame[mainRow,1]<-as.vector(fsr[currRow,"libName"])
      geneFrame[mainRow,2]<-as.vector(fsr[currRow,"term"])
      geneFrame[mainRow,3]<-fsr[currRow,"adjPval"]
      geneFrame[mainRow,4]<-fsr[currRow,"score"]
      geneFrame[mainRow,5]<-genes[gene]
      mainRow <- mainRow + 1
    }
    currRow <- currRow + 1
  }
  return(geneFrame)
}

# libFile needs to have a blank last line.

runEnrichr <- function(libFile, geneFile){
  outputSigFile <- paste(gsub(".csv","",geneFile),"_SIG_EnrichR_Output.csv", sep="")
  outputSigUnpivotFile<- paste(gsub(".csv","",geneFile),"_SIG_EnrichR_UNPIVOT.csv", sep="")
  outputTopFile <- paste(gsub(".csv","",geneFile),"_TOP10_EnrichR_Output.csv", sep="")
  outputTopUnpivot<- paste(gsub(".csv","",geneFile),"_Top10_EnrichR_UNPIVOT.csv", sep="")
  
  geneSet <- fread(geneFile, header=FALSE)
  libSet <- fread(libFile, header=FALSE)
  genes<-as.vector(unlist(geneSet[,1], use.names=F) )
  libs<-as.vector(unlist(libSet[,1], use.names=F))
  
  enr<-enrichr(genes, libs)
  
  fsr <- getSigTerms(enr)
  
  if(dim(fsr)[1] > 0) {
    fwrite(fsr, outputSigFile)
  } else {
    fwrite("No Significant Enrichments", outputSigFile)
    fwrite("No Significant Enrichments", outputTopFile)
  }
}
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS021/Monocle3_gene_modules/')
geneFiles <- Sys.glob('*.txt')
libFile <- "/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/enrichr_libraries_KEGG2019_GObio2018_ARCHS4.txt"
#libFile <- '/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/enrichr_KEGG_GO_ARCHS4.txt'
libSet <- fread(libFile, header=FALSE)
libs<-as.vector(unlist(libSet[,1], use.names=F))
for (geneFile in geneFiles){
  runEnrichr(libFile, geneFile)
}

setwd('/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/Seurat/Aggr4_nocycle_DE_LA-L-A/LA_vs_L')
geneFiles <- Sys.glob('*up*txt')
for (geneFile in geneFiles){
  runEnrichr(libFile, geneFile)
}



# REVIGO can filter enrichr output by collapsing redundant terms----
# http://revigo.irb.hr/


results <- fread('Enrichr/H9noc_D28.vs.GTEx_Cortex.UP.enrichr.GOBioProcess.txt')
results <- results[order(-`Combined Score`)]

#results.go <- results[grepl('GO', libName),]
results.go<- results
results.go[,GOID := tstrsplit(Term, "\\(GO")[2]]
results.go[,GOID := paste0('GO', GOID)]
results.go[,GOID := gsub("\\)", '', GOID)]
results.go
results.go[,term.short := tstrsplit(Term, "\\(")[1]]
results.go
results.go[,term.short := trimws(term.short, which="right")]
results.go
#results.go.clean <- results.go[,c('GOID','adjPval', 'term.short')]
results.go.clean <- results.go[,c('GOID','Adjusted P-value', 'term.short')]

results.go.clean
names(results.go.clean)[2] <- 'adjPval'
#fwrite(results.go.clean, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB008/DE/Partek/CEPT_vs_Y27_Thaw_UP,_SIG_EnrichR_Output_tidy.csv')
fwrite(results.go.clean, 'Enrichr/H9noc_D28.vs.GTEx_Cortex.UP.enrichr.GOBioProcess_output_tidy.csv')



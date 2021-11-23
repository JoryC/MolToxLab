### Libraries
library(tidyverse)
library(edgeR)
library(preprocessCore) # for quantile normalization
library(PMCMRplus)  # for Williams Trend Test

### Folders and file info
setwd("C:\\Users\\nguyenty\\Documents\\R for TDRM\\raw_counts\\")
dataFolder<-"DHEC_nodmso\\"
fNames<-list.files(dataFolder)
sampleNames<-sub("-GQReadsPerGene.out.tab", "", fNames)
sampleNames  # makes sure the doses below match the sampleNames object

# indicate doses
doses<-c(
  rep(0.3,3),
  rep(0.03,3),
  rep(0.003,3),
  rep(0.0003,2),
  rep(3,3)
)

### Import and combine data
rawData<-list()
for(i in 1:length(fNames)){
  rawData[[i]]<-read.table(paste0(dataFolder,fNames[i]),
                           header = FALSE,
                           #header = TRUE,
                           stringsAsFactors = FALSE,
                           sep="\t",
                           strip.white=TRUE
                           )[-c(1:4),]
                           #)
}

countData<-rawData[[1]][,c(1,4)]
for(i in 2:length(rawData)){
  countData<-full_join(countData,
                      rawData[[i]][,c("V1", "V4")], by="V1")
}
colnames(countData)[2:ncol(countData)]<-sampleNames
rownames(countData)<-countData$V1
countData<-countData[,-1]


### Filter low count Genes and Log Transform

# filter out genes with a high poportion of low counts
nonLowProp<-0.7
nonLowRows<-apply(countData, 1, function(x){ (sum(x>10)/length(x))>nonLowProp })
filterData<-countData[nonLowRows,]

# log transform
logFilterData<-log2(filterData +1)


### Quantile normalization

# quantile normalize
quantData<-normalize.quantiles(as.matrix(logFilterData))
rownames(quantData) <- rownames(logFilterData)
colnames(quantData) <- colnames(logFilterData)

# export normalzied data in fastbmd format
saveName<-"no_log_quant_fastbmd_input_DHECv2_nodmso.txt"
s_names<-matrix(sampleNames, nrow=1, dimnames=list(c("#NAME"), NULL))
d_names<-matrix(doses, nrow=1, dimnames=list(c("#CLASS:dose"), NULL))
write.table(s_names, saveName, sep="\t", row.names = TRUE, col.names = FALSE)
write.table(d_names, saveName, sep="\t", row.names = TRUE, col.names = FALSE, append=TRUE)
write.table(quantData, saveName, sep="\t", row.names = TRUE, col.names = FALSE, append=TRUE)


### Filter by Williams Trend Results (import results from BMDExpress)

# import williams trend results from BMDExpress
willResults<-rownames(read.table("DHECv2_nodmso_williams.txt", sep="\t", skip=16, header=TRUE))

# filter quantile normalized data by williams results
wtData<-quantData[rownames(quantData)%in%willResults,]

# export williams filtered data in fastbmd format
saveName<-"wt_fastbmd_input_DHECv2_nodmso.txt"
s_names<-matrix(sampleNames, nrow=1, dimnames=list(c("#NAME"), NULL))
d_names<-matrix(doses, nrow=1, dimnames=list(c("#CLASS:dose"), NULL))
write.table(s_names, saveName, sep="\t", row.names = TRUE, col.names = FALSE)
write.table(d_names, saveName, sep="\t", row.names = TRUE, col.names = FALSE, append=TRUE)
write.table(wtData, saveName, sep="\t", row.names = TRUE, col.names = FALSE, append=TRUE)






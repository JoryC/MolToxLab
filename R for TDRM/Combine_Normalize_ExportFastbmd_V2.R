### Libraries
library(preprocessCore)
library(tidyverse)
library(edgeR)

### Folders and file names
setwd("C:\\Users\\nguyenty\\Documents\\R for TDRM\\raw_counts\\")
dataFolder<-"SANT_nodmso\\"
fNames<-list.files(dataFolder)
sampleNames<-sub("-GQReadsPerGene.out.tab", "", fNames)


### Import and combine data
rawData<-list()
for(i in 1:length(fNames)){
  rawData[[i]]<-read.table(paste0(dataFolder,fNames[i]),
                           header = FALSE,
                           stringsAsFactors = FALSE,
                           sep="\t",
                           strip.white=TRUE
                           )[-c(1:4),]
}

countData<-rawData[[1]][,c(1,4)]
for(i in 2:length(rawData)){
  countData<-full_join(countData,
                       rawData[[i]][,c("V1", "V4")],
                       by="V1")
}

colnames(countData)[2:ncol(countData)]<-sampleNames
rownames(countData)<-countData$V1
countData<-countData[,-1]

# export count data
saveName<-"fastbmd_input_DBP_noDMSO.txt"
write.table(countData, saveName, sep="\t", row.names = TRUE, col.names = TRUE, append=FALSE)


### Quantile normalization


# filter out genes with a high poportion of low counts
nonLowProp<-0.7
nonLowRows<-apply(countData, 1, function(x){ (sum(x>10)/length(x))>nonLowProp })
filterData<-countData[nonLowRows,]

# quantile normalize
logFilterData<-log2(filterData +1)
quantData<-normalize.quantiles(as.matrix(logFilterData))
rownames(quantData) <- rownames(logFilterData)
colnames(quantData) <- colnames(logFilterData)



# export log-transformed and normalzied counts per million data
saveName<-"quant_fastbmd_input_SANT_noDMSO.txt"
write.table(quantData, saveName, sep="\t", row.names = TRUE, col.names = TRUE)
#write.table(logFilterData, saveName, sep="\t", row.names = TRUE, col.names = TRUE)


#Import William's Trend Data
willResults<-rownames(read.table(blank.txt, sep=))
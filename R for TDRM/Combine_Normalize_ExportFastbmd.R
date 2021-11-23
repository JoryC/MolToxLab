### Libraries
library(tidyverse)
library(edgeR)

### Folders and file names
setwd("C:\\Users\\nguyenty\\Documents\\R for TDRM\\raw_counts\\")
folder<-"SANT\\"
fNames<-list.files(folder)
gNames<-sub("-GQReadsPerGene.out.tab", "", fNames)


### Import and combine data
rawData<-list()
for(i in 1:length(fNames)){
  rawData[[i]]<-read.table(paste0(folder,fNames[i]),
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
colnames(countData)[2:ncol(countData)]<-gNames
rownames(countData)<-countData$V1
countData<-countData[,-1]
#write.table(countData, "fastbmd_input_DHEP.txt", sep="\t", row.names=TRUE, col.names=TRUE)


### EdgeR normalization

# set dose groups (1 = control, 2 = low, etc)

gNames  # makes sure the doses below match the gNames object
doses<-c(

  rep(0,4),
  rep(30000,3),
  rep(300,3),
  rep(30,2),
  rep(3,3),
  rep(30000,3)

)

doseFactor<-as.factor(doses)
y<-DGEList(counts<-countData, group=doseFactor)

#Remove low count genes#
keep <- rowSums(cpm(y)>2) >= 3
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples

#Calculate the normalization factor for each column#
yN <- calcNormFactors(y)
yN$samples
plotMDS(yN)

#output normalized counts per million
logcpm <- cpm(y, prior.count=2, log=TRUE)

# export log-transformed and normalzied counts per million data
saveName<-"norm_fastbmd_input_SANT.txt"
dose_names<-matrix(doses, nrow=1)
colnames(dose_names)<-gNames
write.table(dose_names, saveName, sep="\t", row.names = FALSE, col.names = TRUE)
write.table(logcpm, saveName, sep="\t", row.names = TRUE, col.names = FALSE, append=TRUE)

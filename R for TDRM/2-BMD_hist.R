
### LIBRARIES
library(tidyverse)

### FOLDERS, FILES AND FUNCTIONS
setwd("C:\\Users\\nguyenty\\Documents\\R for TDRM")
dataFolder <-"raw_counts\\"
fName<-"all_bmds.txt"

# imprt data
histData<-read.table(paste0(dataFolder,fName), stringsAsFactors = FALSE, header=TRUE, sep="\t")

# load custom mode-detection function
funcFile<-"mode_antimode.R"
source(funcFile)

### VARIABLES AND OPTIONS
histNames<-names(histData)
logBase<-10
min_dense<-0.10  #minimum probability density to be considered a "mode"
bwFun<-"SJ" # the "bandwidth" function to use to determine modes. I've selected the Sheather & Jones (1991) method. see: https://www.ncbi.nlm.nih.gov/pubmed/24885339
xth_gene<-100

# experimental information
lowDose<-c(
  0.0003, #DBC
  0.0003, #DBP
  0.0003, #DHEC
  0.0003, #DHEP
  0.00003, #P1400
  0.0003 #SANT
)

highDose<-c(
  3, #DBC
  3, #DBP
  3, #DHEC
  3, #DHEP
  3, #P1400
  3 #SANT
)

loec<-c(
  10, #DBC
  5, #DBP
  100, #DHEC
  1, #DHEP
  1, #P1400
  5 #SANT 
)

ec50<-c(
  46.2, #DBC
  10.2, #DBP
  370.4, #DHEC
  1000, #DHEP
  22.6, #P1400
  5.4 #SANT 
)


### LISTS TO STORE RESULTS
plotData<-list()
dataModes<-list()
bootModes<-list()
bws<-list()
breaks<-list()

# load data into a list
for(i in 1:length(histNames)){
  plotData[[i]]<-histData[[histNames[i]]]
  plotData[[i]]<-plotData[[i]][!is.na(plotData[[i]])]
  plotData[[i]]<-log(plotData[[i]],logBase)
}


### CALCULATIONS

#calculate bandwidth, breakpoints, modes and antimodes
for(i in 1:length(plotData)){
  bws[[i]]<-bw.SJ(plotData[[i]])
  breaks[[i]]<-seq(from=min(plotData[[i]])-bws[[i]], to=max(plotData[[i]])+bws[[i]], by=bws[[i]])
  dataModes[[i]]<-mode.antimode(plotData[[i]],min.size=min_dense, bw=bwFun)
}

# pull out the first mode
firstMode<-sapply(dataModes, function(x){x$modes[1]})


### PLOTS HISTOGRAMS WITH MODES, DOSES, ETC...

#plot.order<-order(firstMode)
plot.order<-1:6

# auto multi-plot layout
par(mfrow=c(ceiling(length(plot.order)),1)) 

# fixed layout
#par(mfrow=c(6,1))

# plot all graphs
for(i in plot.order){
  #hist(plotData[[i]], breaks=breaks[[i]],  prob=TRUE,    # using algorithm break points and probability density
  hist(plotData[[i]], breaks=50,  prob=FALSE,             # using fixed break points and frequency
       #yaxt="n",
       ylim=c(0,90), # for fixed frequency
       xlim=c(-6,4),
       xlab="log(Dose)",
       main=NULL)
  lines(density(plotData[[i]], bw=bws[[i]]), col=hsv(0.4,1,0.8,0.4), lwd=3)  # pro density curve
  abline(v=log(lowDose[i],logBase), col="blue", lty=3, lwd=3)                # low dose
  abline(v=log(highDose[i],logBase), col="blue", lwd=3)                      # high dose
  abline(v=log(loec[i],logBase), col="orange", lwd=3)                        # loec
  abline(v=log(ec50[i],logBase), col="red", lwd=2)                           # ec50
  #abline(v=dataModes[[i]]$modes[[1]], col="purple", lwd=3)                  # 1st mode
  #abline(v=plotData[[i]][order(plotData[[i]], decreasing = FALSE)][xth_gene], col="orange", lwd=5)  # xth gene
  title(histNames[i], line=-1.5, adj=0.05, bg="white")
}


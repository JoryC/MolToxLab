library(dendextend)

### NOTE: uses "alLData" object from Filter and Normalize script

#covert to cpm
libsize <- allData %>%
  column_to_rownames("sample") %>%
  select(-dose, -chemical) %>%
  apply(1, sum)

pre_cpm <- allData %>%
  column_to_rownames("sample") %>%
  select(-dose, -chemical)

# CPM takes quite a long time....almost 20 minutes! Not any better when "apply" instead of for loop
cpm2<-pre_cpm
for(i in 1:nrow(cpm2)){
  cpm2[i,]<-log2((10^6)*(pre_cpm[i,] + 0.5)/(libsize[i] + 1))
}

#correlation
dataCor <- 1 - cor(t(cpm2), method = "spearman")

#cluster
clustCor <- dataCor %>%
  as.dist() %>%
  hclust(method="average") %>%
  as.dendrogram()

#plot (need to make bottom margin larger so car read labels)
par(mar=c(9,4,1,3))
plot(clustCor)
abline(a=0.1, b=0, col="red", lwd=2) # indicate 0.1 cutoff

library(tidyverse)

fileNames <- list.files("TempRawFiles/Input/")
baseNames <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(fileNames))
loadRaw <- list()

for(i in 1:length(fileNames)){
  loadRaw[[paste0(baseNames[i])]] <- read.delim(file = paste0("TempRawFiles/Input/", fileNames[i]),
                                                header = T)
}

genenames <- read.csv("TempRawFiles/GeneNames/genenames.csv",
                      header = F)

bindgenenames <- list()

for(i in baseNames){
  bindgenenames[[i]] <- cbind(genenames, loadRaw[i]) %>%
    rename("gene" = "V1", "count" = colnames(loadRaw[[i]]))
}

for(i in baseNames){
  write.csv(bindgenenames[[i]], 
            paste0("TempRawFiles/Output/", i, ".csv"),
            row.names = F)
}

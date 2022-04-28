

raw_importer <- function(file_path, file_type, start_phrase, keep_cols = NULL){
  
  #read file names
  filenames <- list.files(file_path, pattern = file_type)
  
  #collect the starting row for each file
  row_start <- vector()
  for(i in 1:length(filenames)){
    row_start[i] <- readLines(paste0(file_path, filenames[i])) %>%
      grepl(start_phrase, .) %>%
      which()
  }
  
  #importing the data
  raw_data <- list()
  for (i in 1:length(filenames)) {
    raw_data[[filenames[i]]] <- suppressWarnings(suppressMessages(  # suppresses annoying parsing error warnings/messages
      read_delim(paste0(file_path, filenames[i]),
                 skip=row_start[i]-1,
                 delim="\t",
                 col_select = keep_cols,
                 show_col_types = FALSE) %>%
        as.data.frame() %>%
        rename_all(make.names)
        
    ))
  }
  return(raw_data)
}


source("BMDExpressFunctions.R")




BMDtest <- raw_importer(file_path = "BMDExpressData/BMD/", 
                 file_type = ".txt", 
                 start_phrase = "Probe Id") %>% 
  lapply(cleanupcolumns)

REACtest <- raw_importer(file_path = "BMDExpressData/REACTOME/",
                         file_type = ".txt",
                         start_phrase = "GO/Pathway/") %>%
  lapply()

GOtest <- raw_importer(file_path = "BMDExpressData/GO_TERM/",
                       file_type = ".txt",
                       start_phrase = "GO/Pathway/")


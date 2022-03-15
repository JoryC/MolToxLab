#### BMD Functions ####
# Clean up columns during import
cleanupcolumns <- function(x){
  x %>%
    select(Probe.Id,
           BMD,
           BMDL,
           BMDU,
           fitPValue)
}

# BMD Filtering
# BMD/BMDL > 20, BMDU/BMDL > 40, BMDU/BMD > 20, BMD > lowdose, BMD < highdose
BMDfiltering <- function(x,
                         BMD.div.BMDL = 20,
                         BMDU.div.BMDL = 40,
                         BMDU.div.BMD = 20,
                         lowdose,
                         highdose,
                         fitP = 0.1){
  x %>%
    filter(
      BMD/BMDL < BMD.div.BMDL,
      BMDU/BMDL < BMDU.div.BMDL,
      BMDU/BMD < BMDU.div.BMD,
      BMD > lowdose/10,
      BMD <= highdose,
      fitPValue > fitP
      )
}

#Get list of BMDs based on lowest median values
BMD_list_extraction <- function(x){
  temp <- x %>% 
    filter(logBMD == min(x[,"logBMD"])) %>%
    select(BMD.List)
  temp2 <- strsplit(temp[,1],";") %>% lapply(as.numeric)
}


####REACTOME Functions####
#Clean up columns during import
cleanupcolumns_reactome <- function(x){
  x %>%
    select(GO.Pathway.Gene.Set.Gene.ID,
           GO.Pathway.Gene.Set.Gene.Name,
           All.Genes..Expression.Data.,
           Genes.That.Passed.All.Filters,
           Fisher.s.Exact.Two.Tail,
           Percentage,
           BMD.Mean,
           BMD.Median,
           BMD.List
           )
}

# Reactome Filtering
# Two Tail Test < 0.05, genes that passed all filters > 3
reactome_filtering <- function(x,
                               p = 0.05,
                               min_gene = 3){
  x %>%
    filter(
      p >= Fisher.s.Exact.Two.Tail,
      min_gene <= Genes.That.Passed.All.Filters
    )
}



#### GENERAL FUNCTIONS ####

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

 
####GO Term Functions####
#Clean up columns
cleanupcolumns_goterm <- function(x){
  x %>%
    select(GO.Pathway.Gene.Set.Gene.ID,
           GO.Level,
           GO.Pathway.Gene.Set.Gene.Name,
           All.Genes..Expression.Data.,
           Genes.That.Passed.All.Filters,
           Fisher.s.Exact.Two.Tail,
           Percentage,
           BMD.Mean,
           BMD.Median,
           BMD.List
    )
}

# GO-term Filtering
# Two Tail Test < 0.05, genes that passed all filters > 3, GO Term level > 5
goterm_filtering <- function(x,
                               p = 0.05,
                               min_gene = 3,
                               min_level = 5){
  x %>%
    filter(
      p >= Fisher.s.Exact.Two.Tail,
      min_gene <= Genes.That.Passed.All.Filters,
      min_level <= GO.Level
    )

}




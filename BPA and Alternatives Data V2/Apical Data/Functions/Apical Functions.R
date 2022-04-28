####combinedanova function####
#calculate anova in list
combinedanova <- function(x, type){
  if(type == "lethal"){
    aov(Survival.Rate ~ Dose, data = x) %>% 
      summary()
  } else if (type == "deform"){
    aov(Deformity.Rate ~ Dose, data = x) %>% 
      summary()
  } else {
    print("Select `lethal` or `deform` as the type")
  }
}

####combineddunnett funtion####
#perform dunnnett's test in list
combineddunnett <- function(df, type){
  if(type == "lethal"){
    DunnettTest(x = df$Survival.Rate, g = df$Dose)
  } else if (type == "deform"){
    DunnettTest(x = df$Deformity.Rate, g = df$Dose)
  } else {
    print("Select `lethal` or `deform` as the type")
  }
}



lethal_combineddunnett <- function(df){
  DunnettTest(x = df$Survival.Rate, g = df$Dose)
}

deform_combineddunnett <- function(df){
  DunnettTest(x = df$Deformity.Rate, g = df$Dose)
}

####summarystats function####
#calculate the mean, counts and SD
summarystats <- function(x, grouping = "Dose", values = "Survival.Rate"){
  group_by(x, Dose) %>%
    summarise(
      count = n(),
      Survival_Rate = mean(Survival.Rate, na.rm = TRUE),
      Survival_StdDev = sd(Survival.Rate, na.rm = TRUE),
      Deform_Rate = mean(Deformity.Rate, na.rm = T),
      Deform_StdDev = sd(Deformity.Rate, na.rm = T)
    )
}

####Multiple plot function####
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


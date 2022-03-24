####Behavioural Data Scripts####

####simi_normalize function####
#normalize data
simi_normalize <- function(x){
  x %>%
    mutate(control_median = unname(tapply(x$endpoint_value, x$is_VC, median)["TRUE"])) %>%
    mutate(endpoint_value_norm = (endpoint_value/control_median)*100-100) %>%
    select(-c(endpoint_value, control_median)) %>%
    mutate(dose = substr(embryo_id, 1, 5)) %>%
    select(-embryo_id)
  # mutate(embryo_id = substr(embryo_id, nchar(embryo_id), nchar(embryo_id)))
}

####dose_replacement function####
#turn dose1-x into actual doses based on metadata
dose_replacement <- function(x, Highdose, foldchange = 10){
  ndosegroups <- 1:(length(unique(x$dose))-1)
  for(i in ndosegroups){
    x$dose[x$dose == paste0("Dose", i)] <- Highdose/foldchange^(i-1)
  }
  x$dose[x$dose == "Contr"] <- 0
  x$dose <- as.numeric(x$dose)
  return(x)
}

####summarystats function####
#calculate the mean, counts and SD
summarystats <- function(x, grouping = "dose", values = "endpoint_value_norm"){
  group_by(x, dose) %>%
    summarise(
      count = n(),
      mean = mean(endpoint_value_norm, na.rm = TRUE),
      sd = sd(endpoint_value_norm, na.rm = TRUE)
    )
}

####combinedanova function####
#calculate anova in list
combinedanova <- function(x){
  aov(endpoint_value_norm ~ dose, data = x) %>% 
    summary()
}

####combineddunnett function####
#perform dunnnett's test in list
combineddunnett <- function(df){
  DunnettTest(x = df$endpoint_value_norm, g = df$dose)
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

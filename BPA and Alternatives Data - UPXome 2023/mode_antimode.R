### MODE DECTOR FUCTION

mode.antimode<-function (x, min.size = 0.1, bw="nrd0", min.bw=NULL) 
{
  #checks
  if (missing(x)) 
    stop("The x argument is required.")
  x <- as.vector(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  if (length(unique(x))<=1) #checks if x is a constant
    return(list(modes = NA, mode.dens = NA, size = 1))
  
  #density function
  dens<-density(x,bw=bw)
  if(!is.null(min.bw)){
    if(dens$bw < min.bw){
      dens<-density(x,bw=min.bw)
    }
  }
  
  #save bw
  mode.bw<-dens$bw
  
  #the difference between each sequential y value
  y.diff<-diff(dens$y)      
  
  #points along the desntiy distribution where y is increasing
  incr<-rep(0,length(y.diff))
  incr[which(y.diff>0)]<-1
  
  #identify points where y changes direction
  begin <- 1
  count <- 1
  for (i in 2:length(incr)) {
    if (incr[i] != incr[i - 1]) {
      count <- count + 1
      begin <- c(begin, i)
    }
  }
  begin <- c(begin, length(incr))
  
  #placeholders for results
  size <- modes <- mode.dens <- rep(0, count/2)
  anti.modes<-rep(0,length(size-1))
  
  #sum of all y values
  sumdens <- sum(dens$y)
  
  ###identify modes
  init <- 1
  
  #if density plot begins on a local max (i.e. a mode) 
  if (incr[1] == 0) {
    size[1] <- sum(dens$y[1:begin[2]])/sumdens
    init <- 2
  }
  
  #identfy all modes, size and max density
  j <- init
  for (i in init:length(size)) {
    size[i] <- sum(dens$y[begin[j]:begin[j + 2]])/sumdens
    temp.x<- dens$x[begin[j]:begin[j + 2]]
    temp.y<- dens$y[begin[j]:begin[j + 2]]
    highs<-which(temp.y==max(temp.y))
    if(median(highs)%%1>0){
      high.point<-median(highs[-length(highs)])
    }else{
      high.point<-median(highs)
    }
    modes[i] <- temp.x[high.point]
    mode.dens[i] <- temp.y[high.point]
    j <- j + 2
  }

  #identify and remove modes smaller than min.size
  if (any(size < min.size)) {
    modes <- modes[-which(size < min.size)]
    mode.dens <- mode.dens[-which(size < min.size)]
    size <- size[-which(size < min.size)]
  }
  
  #identify anti-mode between modes
  if(length(modes)>1){
    anti.modes<-rep(0,length(size)-1)
    for(i in 1:length(anti.modes)){
      m1<-which(dens$y==mode.dens[i])
      m1<-m1[m1%in%begin]
      m2<-which(dens$y==mode.dens[i+1])
      m2<-m2[m2%in%begin]
      temp.x <- dens$x[m1:m2]
      temp.y <- dens$y[m1:m2]
      lows<-which(temp.y==min(temp.y))
      if(median(lows)%%1>0){
        low.point<-median(lows[-length(lows)])
      }else{
        low.point<-median(lows)
      }
      anti.modes[[i]]<-temp.x[low.point]
    }
  }else{
    anti.modes<-NULL
  }
    
  return(list(modes = modes, mode.dens = mode.dens, size = size, anti.modes=anti.modes, bw=mode.bw))
}






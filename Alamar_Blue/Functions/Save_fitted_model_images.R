library(png)
library(grid)
library(gridExtra)
library(dplyr)

png_list <- list()
for (i in 1:length(chemicalNames)) {
  p <- readPNG(paste0(
    "Output/",
    chemicalNames[i], "/", chemicalNames[i], 
    "_hill_fitted_", "z_score", "_", dir, "_dir", params$sample_size, "_samples.png"
  ))
  png_list[[chemicalNames[i]]] <- p
  rm(p)
}

yleft <-
  grid::textGrob(paste0("Z-Score"),
                 rot = 90,
                 gp = grid::gpar(fontsize = 24))
bottom <- grid::textGrob("Dose (mg/L)", gp = grid::gpar(fontsize = 24))

png(
  paste0(
    "Output/Images/",
    "all_chems_hill_fitted_", "z_score", "_", params$direction, "_dir_",
    params$sample_size,
    "_samples.png"
  ),
  units = "px",
  width = 1920,
  height = 1080,
  bg = "white"
)

gridExtra::grid.arrange(
  rasterGrob(png_list[[1]]),
  rasterGrob(png_list[[2]]),
  rasterGrob(png_list[[3]]),
  rasterGrob(png_list[[4]]),
  rasterGrob(png_list[[5]]),
  rasterGrob(png_list[[6]]),
  rasterGrob(png_list[[7]]),
  rasterGrob(png_list[[8]]),
  rasterGrob(png_list[[9]]),
  rasterGrob(png_list[[10]]),
  rasterGrob(png_list[[11]]),
  rasterGrob(png_list[[12]]),
  rasterGrob(png_list[[13]]),
  rasterGrob(png_list[[14]]),
  rasterGrob(png_list[[15]]),
  rasterGrob(png_list[[16]]),
  rasterGrob(png_list[[17]]),
  rasterGrob(png_list[[18]]),
  rasterGrob(png_list[[19]]),
  rasterGrob(png_list[[20]]),
  rasterGrob(png_list[[21]]),
  rasterGrob(png_list[[22]]),
  rasterGrob(png_list[[23]]),
  rasterGrob(png_list[[24]]),
  rasterGrob(png_list[[25]]),
  rasterGrob(png_list[[26]]),
  rasterGrob(png_list[[27]]),
  rasterGrob(png_list[[28]]),
  rasterGrob(png_list[[29]]),
  plot(x= 0,y=0, ann=FALSE, axes=FALSE, col="white"),
  left = yleft,
  bottom = bottom,
  nrow = 5,
  ncol = 6
)

dev.off()
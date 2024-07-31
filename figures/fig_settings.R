fsize=9
require(dplyr)
require(matrixTests)
require(ggplot2)
require(ggplotify)
require(patchwork)
require(data.table)
require(parallel)
require(devtools)
require(annmatrix)
require(pROC)
require(lomb)
load_all()
theme_set(theme_classic()) 

show_temp_plt=function(plt,plt_width,plt_height){
  if (interactive()){
    plt_path <- tempfile(fileext = ".png")
    ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
           dpi = 96)
    viewer <- getOption("viewer")
    viewer(plt_path)
  } 
}

pub_qual=F
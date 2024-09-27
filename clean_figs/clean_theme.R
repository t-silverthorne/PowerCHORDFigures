clean_theme=function(){
  list(
    theme_classic(),
    theme(strip.background = element_blank()),
    theme(strip.text = element_text(
      size = 7, color = "black",
      margin = margin(b = 1, t = 1)
    )),
    theme(
    plot.margin = margin(0,0,0,0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)),
    theme(axis.text = element_text(size = 7)),
    theme(text=element_text(size=7))
  )
}
show_temp_plt=function(plt,plt_width,plt_height){
  if (interactive()){
    plt_path <- tempfile(fileext = ".png")
    ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
           dpi = 96)
    viewer <- getOption("viewer")
    viewer(plt_path)
  } 
}

require(ggplot2)
require(dplyr)
require(matrixTests)
require(ggplot2)
require(data.table)
require(tidyr)
require(patchwork)
require(parallel)
require(devtools)
require(annmatrix)
devtools::load_all()
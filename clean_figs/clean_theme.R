require(matrixTests)
require(matrixTests)
require(ggplot2)
require(data.table)
require(tidyr)
require(patchwork)
require(ggplot2)
require(dplyr)
pub_qual=F
mar=2
clean_theme=function(){
  list(
    theme_classic(),
#    theme(strip.background = element_blank()),
    theme(strip.text = element_text(
      size = 7, color = "black",
      margin = margin(b = mar,r =mar,t=mar,l=mar)
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

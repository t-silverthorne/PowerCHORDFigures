# always need
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

# might need
plt = plt + labs(x=element_text('x variable'),
                 y=element_text('y variable'),
                 color='z variable')


# always need
plt=plt+theme(text=element_text(size=fsize))
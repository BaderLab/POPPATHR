theme_Publication <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
            #   axis.title.x=element_blank(),
            #   axis.text.x = element_text(angle=45, hjust=1),
               axis.text = element_text(),
            #   axis.text.x = element_blank(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
            #   axis.ticks.x = element_blank(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
            #   legend.position = "right",
               legend.direction = "horizontal",
            #   legend.direction = "vertical",
               legend.key.size = unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
          #     legend.title = element_text(face="italic"),
               legend.title = element_blank(),
               strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold"),
               plot.margin = unit(c(10,5,5,5),"mm"),
            #   plot.margin = unit(c(0.5,4,0.5,8),"cm")
          ))

}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill", "Publication",
          manual_pal(values=c("#fb9a99", "#386cb0", "#984ea3", "#ffff33",  "#fdb462",
                              "#7fc97f", "#ef3b2c", "#662506", "#a6cee3"
                              )), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour", "Publication",
          manual_pal(values=c("#fb9a99","#386cb0", "#984ea3", "#ffff33",  "#fdb462",
                              "#7fc97f", "#ef3b2c", "#662506", "#a6cee3"
                              )), ...)

}

#fdb462 - dark yellow
#984ea3 - purple
#386cb0 - blue
#7fc97f - green
#ef3b2c - red
#a6cee3 - light blue
#ffff33 - bright yellow
#662506 - brown
#fb9a99 - pink

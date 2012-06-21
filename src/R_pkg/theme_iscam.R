# Default plot theme for iscam in ggplot2

theme_iscam <- function (base_size = 12, base_family = "") 
{
    structure(list(
        axis.line         = theme_blank(),
        axis.text.x       = theme_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", vjust = 1, angle = 00),
        axis.text.y       = theme_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1, angle = 00),
        axis.ticks        = theme_segment(colour = "grey50"),
        axis.title.x      = theme_text(size = base_size, colour = "grey50", face="bold"),
        axis.title.y      = theme_text(size = base_size, colour = "grey50", face="bold", angle = 90),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks.margin = unit(0.10, "cm"),

        legend.background = theme_rect(colour= NA), 
        legend.key        = theme_rect(fill = "grey95", colour = "white"),
        legend.key.size   = unit(1.2, "lines"),
        legend.text       = theme_text(size = base_size * 0.7, colour = "grey50"),
        legend.title      = theme_text(size = base_size * 0.7, colour = "grey50", hjust = -1),
        legend.position   = "right",

        panel.background  = theme_rect(fill = "grey90", colour = NA), 
        panel.border      = theme_blank(), 
        panel.grid.major  = theme_line(colour = "white"),
        panel.grid.minor  = theme_line(colour = "grey95", size = 0.25),
        panel.margin      = unit(0.25, "lines"),

        strip.background  = theme_rect(fill = "grey80", colour = NA), 
        strip.label       = function(variable, value) value, 
        strip.text.x      = theme_text(size = base_size * 0.8),
        strip.text.y      = theme_text(size = base_size * 0.8, angle = -90),

        plot.background   = theme_rect(colour = NA),
        plot.title        = theme_text(size = base_size * 1.0, colour = "grey50"),
        plot.margin       = unit(c(1, 1, 0.5, 0.5), "lines")
    ), 
    class = "options")
}

# A simple function for creating transparent colors
# Author: Nathan Stephens (hacks package)
colr <- 
function(col.pal=1,a=1)
{
    col.rgb<-col2rgb(col.pal)/255
    rgb(t(col.rgb),alpha=a)
}
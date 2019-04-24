#!/use/bin/env R
require(grid)
require(ggplot2)

### theme definitions ######################################################
thm = list()

thm$commonTheme = theme_grey() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.key       = element_blank(), 
    legend.background= element_rect(fill="white", size=0), 
    panel.background = element_blank(), 
    panel.border     = element_blank(), 
    strip.background = element_blank())

thm$pdfHeight = 6
thm$pdfWidth  = 8
thm$pdfTheme  = theme(
        axis.line       = element_line(size = 0.7, color = "black"),
        axis.ticks      = element_line(color = "darkgrey",size=0.7),
        axis.ticks.length = unit(0.15,"cm"),
        title           = element_text(size=16, face="bold"),
        text            = element_text(size = 10),
        axis.text       = element_text(size = 12), 
        axis.title      = element_text(face="plain"), 
        legend.title    = element_text(size = 12, face = "bold"), 
        legend.text     = element_text(size = 12),
        panel.spacing   = unit(1,"cm"))

thm$pngWidth   = 1600
thm$pngHeight  = 1200
thm$pngTheme   = theme(
        axis.line       = element_line(size = 2, color = "black"),
        axis.ticks      = element_line(color="darkgrey",size=2),
        axis.ticks.length = unit(0.4,"cm"),
        title           = element_text(size=44, face="bold"),
        text                  = element_text(size = 30),
        axis.text             = element_text(size = 36), 
        axis.title            = element_text(face="plain"), 
        legend.title          = element_text(size = 36, face = "bold"), 
        legend.text           = element_text(size = 36),
        legend.key.height     = unit(3,"line"),
        legend.key.width     = unit(3,"line"),
        plot.margin    = unit(c(1.6,1.6,1.6,1.6),"cm"))

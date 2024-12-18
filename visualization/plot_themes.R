add_fonts <- function() {
  
  sysfonts::font_add(family = "Avenir Next",
                     regular = "fonts/AvenirNext-Regular.ttf",
                     bold = "fonts/AvenirNext-DemiBold.ttf")
  sysfonts::font_add(family = "Georgia",
                     regular = "fonts/Georgia.ttf")
   showtext::showtext_auto()
  
  cowplot::set_null_device("png")
  
}

theme_manuscript2 <- function() {
  theme_bw() %+replace%
    theme(
      text = element_text(family = "Avenir Next"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7),
      axis.ticks = element_line(color = "#BDBDBD"),
      axis.line = element_line(color = "#BDBDBD"),
      legend.background = element_rect(fill = "transparent",
                                       color = "transparent"),
      legend.margin = margin(0, 0, 0, 0, "pt"),
      legend.spacing.x = unit(8, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.key = element_rect(fill = "transparent",
                                color = "transparent"),
      
      legend.key.width = unit(20, "pt"),
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.box.margin = margin(0, 0, 0, 0, "pt"),
      legend.box.background = element_rect(fill = "transparent",
                                           color = "transparent"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
      plot.title = element_text(size = 9),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      plot.margin = margin(t = 13, r = 10, b = 5.5, l = 10, "pt"),
      complete = TRUE
    )
}

theme_manuscript <- function() {
  theme_manuscript2() %+replace%
    theme(legend.text = ggtext::element_markdown(size = 9))
}


theme_presentation <- function() {
  theme_manuscript() %+replace%
    theme(text = element_text(family = "Avenir Next", size = 18),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.background = element_rect(fill = "#F2F2F2", color = "#F2F2F2"),
          plot.title = element_text(size = 18, hjust = 0.5),
          complete = TRUE)
}

theme_poster <- function() {
  theme_bw() %+replace%
    theme(
      text = element_text(family = "Avenir Next"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = NA, fill = NA),
      plot.background = element_rect(fill = "#F0F0F0", color = "#F0F0F0"),
      plot.title = element_text(size = 20, hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5),
      plot.margin = margin(t = 25, r = 10, b = 5.5, l = 10, "pt"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 20),
      axis.line = element_line(color = "#000000"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.background = element_rect(fill = "#F0F0F0", color = "#F0F0F0"),
      legend.key = element_rect(fill = "#F0F0F0", color = "#F0F0F0"),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 18),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.margin = margin(0, unit = "mm"),
      complete = TRUE
    )
}

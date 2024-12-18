get_legend <- function(plt, position = "bottom") {
  component <- paste0("guide-box-", position)
  lh <- cowplot::get_plot_component(plt, component)
  return(lh)
}

add_tag <- function(label) {
  ggplot2::geom_text(data = tibble(x = 0, y = 1, label = label),
            aes(x = x, y = y, label = label),
            color = "black", family = "Avenir Next", fontface = "bold",
            size = (10 / .pt), hjust = -1, vjust = 1.2)
}

add_title <- function(label, size = 10) {
  ggplot2::geom_text(data = tibble(x = 0.5, y = 1, label = label),
                     aes(x = x, y = y, label = label),
                     color = "black", family = "Avenir Next", fontface = "bold",
                     size = (size / .pt), hjust = 0.4, vjust = 1.2)
}

label_generations <- function(generation) {
  return(paste0("t = ", generation))
}

define_layout <- function(subplots, nrow = NULL, ncol = NULL) {
  
  n_subplots <- length(subplots)
  
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(n_subplots))
    nrow <- ceiling(n_subplots/ncol)
  }
  
  if (is.null(ncol)) { ncol <- ceiling(n_subplots/nrow) }
  if (is.null(nrow)) { nrow <- ceiling(n_subplots/ncol) }
  
  return(c(nrow, ncol))
  
}

add_legend <- function(plot, lh = NULL, bg = "#FFFFFF") {
  if (!is.null(lh)) {
    lh_height <- sum(lh$heights)
    plot_height <- unit(1, "npc") - lh_height
    plot <- gridExtra::arrangeGrob(grobs = list(plot, lh),
                                   heights = grid::unit.c(plot_height,
                                                          lh_height))
    plot <- cowplot::plot_grid(plot) +
      theme(plot.background = element_rect(fill = bg, color = bg))
  }
  return(plot)
}

add_colorbar <- function(plot, cb = NULL, bg = "#FFFFFF") {
  
  if (!is.null(cb)) {
    cb_width <- sum(cb$widths)
    plot_width <- unit(1, "npc") - cb_width
    plot <- gridExtra::arrangeGrob(grobs = list(plot, cb),
                                   widths = grid::unit.c(plot_width,
                                                          cb_width))
    plot <- cowplot::plot_grid(plot) +
      theme(plot.background = element_rect(fill = bg, color = bg))
  }
  
  return(plot)
}

assemble_panel <- function(plots, nrow = NULL, ncol = NULL, align_y = FALSE,
                           lh = NULL, cb = NULL, title = NULL, bg = "#FFFFFF") {
  
  if (is.ggplot(plots)) {
    plots <- list(plots)
  }
  
  dims <- define_layout(plots, nrow, ncol)
  nrow <- dims[1]
  ncol <- dims[2]
  
  y <- sapply(plots, function(x) { return(layer_scales(x)$y$range$range)})
  y <- c(min(y), max(y))
  
  for (p in seq_along(plots)) {
    
    m <- plots[[p]]$theme$plot.margin
    if ((p + 1) %% ncol == 1) { m[2] <- unit(5.5, "pt") }
    if (p %% ncol == 1) { m[4] <- unit(5.5, "pt") }
    plots[[p]]$theme$plot.margin <- m
    
    if (align_y) {
      plots[[p]] <- plots[[p]] +
        coord_cartesian(ylim = y)
      }
  }
  
  panel <- cowplot::plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)
  
  if (!is.null(title)) {
    panel <- panel + add_title(title)
  }
  panel <- add_legend(panel, lh = lh, bg = bg)
  panel <- add_colorbar(panel, cb = cb, bg = bg)

  return(panel)
}

assemble_figure <- function(panels, nrow = NULL, ncol = NULL, lh = NULL,
                            cb = NULL, add_tags = TRUE, bg = "#FFFFFF", ...) {
  
  dims <- define_layout(panels, nrow, ncol)
  nrow <- dims[1]
  ncol <- dims[2]
  
  col <- seq_along(panels) %% ncol
  col[col == 0] <- ncol
  row <- ceiling(seq_along(panels) / ncol)
  
  for (p in seq_along(panels)) {
    panels[[p]] <- cowplot::plot_grid(panels[[p]])
    
    pos_col <- p %% ncol
    pos_row <- ceiling(p / ncol)
    
    t <- case_when(pos_row == 1 ~ 0, .default = 8)
    r <- case_when(pos_col == 0 ~ 0, .default = 8)
    b <- case_when(pos_row == nrow ~ 0, .default = 8)
    l <- case_when(pos_col == 1 ~ 0, .default = 8)
    panels[[p]] <- panels[[p]] +
      theme(plot.margin = margin(t = t, r = r, b = b, l = l, unit="pt"))
  
    if (add_tags) {
      
      panels[[p]] <- panels[[p]] +
        add_tag(letters[p])
    }
  }
  
  figure <- cowplot::plot_grid(plotlist = panels, ncol = ncol, nrow = nrow)
  
  figure <- add_legend(figure, lh = lh, bg = bg)
  figure <- add_colorbar(figure, cb = cb, bg = bg)
  
  return(figure)
}

save_figure <- function(figure_directory, fig = last_plot(),
                        counter, supplementary = FALSE,
                        height = 6.5, width = 6.5) {
  counter <- counter + 1
  figure_name <- case_when(supplementary ~ paste0("FigureS", counter),
                           .default = paste0("Figure", counter))
  filename <- paste0(figure_directory, "/", figure_name, ".pdf")
  cowplot::save_plot(filename, fig, base_height = height, base_width = width)
  return(counter)
}

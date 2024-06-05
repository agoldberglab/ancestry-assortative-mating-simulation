library(tidyverse)
library(ggtext)

source("visualization/load_data.R")
source("visualization/plot_themes.R")
source("visualization/plot_utils.R")
source("visualization/plotting_functions.R")

add_fonts()

update_geom_defaults("path", aes(linewidth = 0.4))
update_geom_defaults("point", aes(size = 0.6))

lh_linewidth <- 0.6
lh_point_size <- 1

counter <- 0
figure_directory <- "figures"

##############################################################################

# Figure 1

models_shown <- c("exponential_decay", "exponential_decay_normalized",
                  "normal_distribution")

x1_shown <- c(0, 0.5)

panels <- list()

for (m in seq_along(models_shown)) {
  
  plts <- list()
  
  for (x in seq_along(x1_shown)) {
    
    plt <- ancestry_psi %>%
      filter(model == models_shown[m], x1 == x1_shown[x]) %>%
      plot_ancestry_psi(., titles, values, labels) +
      theme_manuscript() +
      theme(plot.margin = margin(13, 13.75, 5.5, 13.75, "pt"),
            axis.title.y = element_markdown())
    
    if (is.null(plt$guides$guides$colour)) {
      plt <- plt +
        guides(color =
                 guide_legend(byrow = TRUE,
                              override.aes = list(linewidth = lh_linewidth)))
    }
    
    if (is.null(plt$guides$guides$linetype)) {
      plt <- plt +
        guides(linetype =
                 guide_legend(byrow = TRUE,
                              override.aes = list(linewidth = lh_linewidth)))
    }
    
    lh <- get_legend(plt)
    
    plts[[x]] <- plt +
      theme(legend.position = "none")
    
  }
  
  panels[[m]] <- assemble_panel(plts, ncol = 2, align_y = TRUE, lh = lh,
                                title = labels$model[models_shown[m]])
}

plt <- social_psi %>%
  plot_social_psi(., titles, values, labels) +
  theme_manuscript() +
  theme(axis.text.x = element_text(size = 10),
        plot.margin = margin(22, 5.5, -20, 5.5, "pt"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.margin = margin(-5, 0, 25, 0, "pt"),
        axis.title.y = ggtext::element_markdown())

panels[[4]] <- assemble_panel(plt, title = labels$model['subpop'])

fig <- assemble_figure(panels)
counter <- save_figure(figure_directory, fig, counter, height = 4)

##############################################################################

# Figure 2

models_shown <- models
thetas_shown <- c(seq(2, 10, by = 2), seq(25, 100, by = 25))
generations_shown <- seq(0, 20)
migration_shown <- 0

plts <- list()

for (t in seq_along(thetas_shown)) {
  
  plt <- summary_data %>%
    filter(model %in% models_shown,
           theta %in% c(1, thetas_shown[t]),
           generation %in% generations_shown,
           migration %in% migration_shown) %>%
    plot_variable_over_time(., "parent_corr", "model",
                            titles, values, labels) +
    theme_manuscript() +
    coord_cartesian(ylim = c(0, 1)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                override.aes = list(linewidth = lh_linewidth)),
           linetype = "none") +
    labs(y = "Corr. between mates", 
         title = labels$theta[as.character(thetas_shown[t])]) +
    theme(plot.title = element_markdown())
  
  lh <- get_legend(plt)
  
  plts[[t]] <- plt +
    theme(legend.position = "none")
}
  
panels <- list()

panels[[1]] <- assemble_panel(plts[1:3], ncol = 3)
panels[[2]] <- assemble_panel(plts[4:6], ncol = 3)
panels[[3]] <- assemble_panel(plts[7:9], ncol = 3)

fig <- assemble_figure(panels, nrow = 3, lh = lh, show_tags = FALSE)
counter <- save_figure(figure_directory, fig, counter)

##############################################################################

# Figure 3

models_shown <- c("exponential_decay_normalized", "subpop")

fig <- plot_correlation_figure(example_pedigree, geom = "hex",
                               examples = c(3, 3),
                        models = models_shown,
                        n = c(10000, 100))

counter <- save_figure(figure_directory, fig, counter, height = 5)

##############################################################################

# Figure 4

models_shown <- c("exponential_decay", "exponential_decay_normalized",
                  "normal_distribution", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2))
generations_shown <- seq(0, 20)
migration_shown <- 0

plts <- list()

for (m in seq_along(models_shown)) {
  
  plt <- summary_data %>%
    filter(model %in% c("random", models_shown[m]),
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown) %>%
    plot_variable_over_time(., "ancestry0_var", "theta",
                            titles, values, labels) +
    guides(linetype = "none",
           color = guide_legend(nrow = 1, byrow = TRUE,
                                override.aes = 
                                  list(linewidth = lh_linewidth))) +
    theme_manuscript()
  
  lh <- get_legend(plt)
  
  plt <- plt + theme(legend.position = "none")
  
  plts[[m]] <- assemble_panel(plt, title = labels$model[models_shown[m]])
}

fig <- assemble_figure(plts, lh = lh)
counter <- save_figure(figure_directory, fig, counter, height = 6)

##############################################################################

# Figure 5

models_shown <- c("random", "exponential_decay_normalized", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2))
generations_shown <- 20
migration_shown <- 0

panels <- list()

panels[[1]] <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown) %>%
  plot_correlation_vs_B(., titles, values, labels) +
  theme_manuscript() +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(0, 14, 0, 0, unit = "pt")))
        
plt <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown) %>%
  ggplot(aes(x = B, y = ancestry0_var, color = model)) +
  geom_point() +
  scale_color_manual(values = values$model, labels = labels$model) +
  labs(x = titles['expressed_mating_bias'], y = titles['ancestry0_var'],
       title = "Generation 20") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(size = lh_point_size))) +
  coord_cartesian(ylim = c(0,0.008)) +
  theme_manuscript()

lh <- get_legend(plt)

panels[[2]] <- plt +
  theme(legend.position = "none")

fig <- assemble_figure(panels, ncol = 2, lh = lh)
counter <- save_figure(figure_directory, fig, counter, height = 3.5)

##############################################################################

# Figure 6

models_shown <- c("random", "exponential_decay_normalized", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2))
generations_shown <- 20
migration_shown <- 0

data <- summary_data %>%
  filter(!is.na(example), generation %in% generations_shown) %>%
  arrange(example, model) %>%
  mutate(x = c(0.78, 1.78, 2.23, 2.78, 3.23, 3.78, 4.23),
  lab = paste(symbol_theta, " =", theta))

panels <- list()

xticks <- unique(round(data$parent_corr, digits = 2))
xticks <- sprintf(fmt = "%0.02f", xticks + 0)

panels[[1]] <- data %>%
  ggplot() +
  geom_boxplot(aes(x = example, ymin = tracts_min, lower = tracts_q25,
                   middle = tracts_q50, upper = tracts_q75, ymax = tracts_max,
                   fill = model),
               stat='identity', width = 0.8,
               position = position_dodge(preserve = "single", width = 0.9)) +
  geom_richtext(aes(x = x, label = lab), y = -0.5, family = "Avenir Next", 
            size = 5.5 * 0.36, fill = NA, label.colour = NA) +
  scale_fill_manual(values = values$model, labels = labels$model) +
  scale_x_discrete(labels = paste("r = ", xticks)) +
  theme_manuscript() +
  labs(x = "", y = "Tract length, cM", title = "Generation 20") +
  theme(legend.position = "none")

plt <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown) %>%
  ggplot(mapping = aes(x = parent_corr, y = tracts_q50, color = model,
                       fill = model)) +
  geom_point(key_glyph = draw_key_rect) +
  theme_manuscript() +
  scale_color_manual(values = values$model, labels=labels$model) +
  scale_fill_manual(values = values$model, labels=labels$model) +
  labs(x = titles['parent_corr'], y = titles['tract_median'],
       title = "Generation 20") +
  theme(legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(0.8, "lines"))

lh <- get_legend(plt)

panels[[2]] <- plt + theme(legend.position = "none")

fig <- assemble_figure(panels, lh = lh)
counter <- save_figure(figure_directory, fig, counter, height = 3.5)

##############################################################################

# Figure 7

models_shown <- c("random", "exponential_decay_normalized", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2))
generations_shown <- 20
migration_shown <- 0

data <- summary_data %>%
  filter(!is.na(example), generation %in% generations_shown) %>%
  arrange(example, model) %>%
  mutate(x = c(0.78, 1.78, 2.23, 2.78, 3.23, 3.78, 4.23),
         lab = paste(symbol_theta, " =", theta)) %>%
  select(x, tracts_mean, tracts_exp_mean, tracts_q95, tracts_exp_q95,
         model, theta, parent_corr, lab) %>%
  pivot_longer(cols = tracts_mean:tracts_exp_q95, values_to = "length") %>%
  mutate(stat = case_when(name %in% c("tracts_mean", "tracts_exp_mean") ~ "mean",
                          .default = "q95"),
         group = case_when(name %in% c("tracts_mean", "tracts_q95") ~ "observed",
                           .default = "expected")) %>%
  mutate(group = factor(group, levels = c("observed", "expected")))

values$tract_length <- c("observed" = "solid", "expected" = "11")
labels$tract_length <- c("observed" = "Observed",
                         "expected" = "Expected under exponential fit")

xticks <- unique(round(data$parent_corr, digits = 2))
xticks <- sprintf(fmt = "%0.02f", xticks + 0)

stats <- c("mean", "q95")
ylabs <- c("Mean tract length, cM", "95th percentile length, cM")

lambda <- 21 * 0.5
expected <- c(1/lambda*100, -log(1-0.95)/lambda*100)

panels <- list()

for (p in seq(2)) {
  plt <- data %>%
    filter(stat == stats[p]) %>%
    ggplot(aes(x = x - 0.1, xend = x + 0.1, y = length, color = model,
               linetype = group)) +
    geom_segment() +
    geom_hline(aes(yintercept = expected[p]), linetype = "dashed",
               color = "#737373", linewidth = 0.2, alpha = 0.7) +
    geom_richtext(aes(x = x, label = lab), y = -0.5, family = "Avenir Next", 
                  size = 5.5 * 0.36, colour="black",
                  fill = NA, label.colour = NA) +
    scale_color_manual(values = values$model, labels = labels$model) +
    scale_linetype_manual(values = values$tract_length,
                          labels = labels$tract_length) +
    labs(x = "", y = "Mean tract length, cM") +
    coord_cartesian(ylim = c(0, 45)) +
    guides(color = "none",
           linetype = guide_legend(nrow = 1, byrow = TRUE,
                                   override.aes =
                                     list(linewidth = lh_linewidth))) +
    theme_manuscript()
  
  lh <- get_legend(plt)
  
  panels[[p]] <- plt + theme(legend.position = "none")
  
}

panelA <- assemble_figure(panels[1:2], lh = lh) 

examples_shown <- c(2, 3, 4)
for (e in seq_along(examples_shown)) {
  plt <- plot_estimated_vs_true_time(examples_shown[e])
  
  panels[[e+2]] <- plt +
    theme(legend.position = "none")
}

for (p in seq_along(panels)) {
  panels[[p]] <- assemble_panel(panels[[p]]) +
    add_tag(LETTERS[p])
}

lh <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown) %>%
  ggplot(mapping = aes(x = parent_corr, y = tracts_q50, color = model,
                       fill = model)) +
  geom_point(key_glyph = draw_key_rect) +
  scale_color_manual(values = values$model, labels=labels$model) +
  scale_fill_manual(values = values$model, labels=labels$model) +
  theme_manuscript() +
  theme(legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(0.8, "lines"))

lh <- get_legend(lh)

panelB <- assemble_figure(panels[3:5], lh = lh, nrow = 1, show_tags = FALSE)

fig <- assemble_figure(list(panelA, panelB), ncol = 1, show_tags = FALSE)
counter <- save_figure(figure_directory, fig, counter, height = 5)
##############################################################################

# Figure 8

models_shown <- c("exponential_decay_normalized", "subpop")
thetas_shown <- c(1, 2, 4, 10)
generations_shown <- seq(0, 20)
migration_shown <- c(0, 0.01)

panels <- list()

for (m in seq_along(models_shown)) {
  
  plt <- summary_data %>%
    filter(model %in% c("random", models_shown[m]),
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown) %>%
    plot_variable_over_time(., "parent_corr", "theta",
                            titles, values, labels) +
    guides(color = guide_legend(byrow = TRUE,
                                override.aes =
                                  list(linewidth = lh_linewidth)),
           linetype = guide_legend(byrow = TRUE,
                                   override.aes =
                                     list(linewidth = lh_linewidth))) +
    theme_manuscript() +
    theme(legend.box = "vertical")
  
  lh <- get_legend(plt)
  
  plt <- plt + theme(legend.position = "none")
  
  panels[[m]] <- assemble_panel(plt) +
    add_title(labels$model[models_shown[m]])

}

migration_shown <- 0.01

for (m in seq_along(models_shown)) {
  
  plt <- summary_data %>%
    filter(model %in% c("random", models_shown[m]),
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown) %>%
    plot_variable_over_time(., "p_migrant_pairs", "theta",
                            titles, values, labels) +
    coord_cartesian(ylim = c(0, 0.41)) +
    theme_manuscript() +
    theme(legend.position = "none")
  
  panels[[m+2]] <- assemble_panel(plt) +
    add_title(labels$model[models_shown[m]])

}

fig <- assemble_figure(panels, lh = lh, ncol = 2)
counter <- save_figure(figure_directory, fig, counter, height = 5.5)

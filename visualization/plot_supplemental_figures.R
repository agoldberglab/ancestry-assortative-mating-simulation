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
# Figure S1

simulations_shown <- unique(example_psi_ratio$simulation)

generations_shown <- seq(5, 20, by = 5)
colors <- set_names(c("#FBB4B9", "#F768A1", "#C51B8A", "#7A0177"),
                    generations_shown)

panels <- list()

for (s in seq_along(simulations_shown)) {
  
  data <- example_psi_ratio %>%
    filter(simulation == simulations_shown[s],
           generation %in% generations_shown)
  
  plt <- data %>%
    plot_psi_over_time(., generations_shown, colors, labels) +
    guides(color = guide_legend(override.aes = list(size = lh_point_size))) +
    theme_manuscript() +
    theme(axis.title.y = element_markdown())
  
  lh <- get_legend(plt)
  
  plt <- plt + theme(legend.position = "none")
  
  title <- labels$short_model[unique(as.character(data$model))]
  panels[[s]] <- assemble_panel(plt, title = title)
  
}

fig <- assemble_figure(panels, lh = lh, ncol = 3)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 2.5)

##############################################################################

# Figure S2 / Figure S3 / Figure S4 / Figure S5

models_shown <- c("exponential_decay", "exponential_decay_normalized",
                  "normal_distribution", "subpop")
thetas_shown <- 5
generations_shown <- seq(1, 20)
migration_shown <- 0

simulations_shown <- binned_ancestry0 %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown,
         seed == 2547432215832) %>%
  pull(simulation) %>%
  unique()

for (s in simulations_shown) {
  fig <- plot_ancestry_histogram_over_time(s, generations_shown)
  counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                         height = 6.5)
}

#############################################################################

# Figure S6

models_shown <- c("exponential_decay", "exponential_decay_normalized",
                  "normal_distribution", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2), seq(25, 100, by = 25))
generations_shown <- seq(0, 20)
migration_shown <- 0

plts <- list()
  
for (m in seq_along(models_shown)) {
  
  data <- summary_data %>%
    filter(model %in% c("random", models_shown[m]),
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown)
  
  plt <- plot_variable_over_time(data, "parent_corr", "theta", titles,
                                 values, labels) +
    guides(color = guide_legend(byrow = TRUE,
                                override.aes = list(linewidth = lh_linewidth)),
           linetype = "none") +
    theme_manuscript()
  
  lh <- get_legend(plt)
  plt <- plt +
    theme(legend.position = "none")
  
  plts[[m]] <- assemble_panel(plt, title = labels$model[[models_shown[m]]])
  
}

fig <- assemble_figure(plts, lh = lh)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE)

#############################################################

# Figure S7

models_shown <- c("exponential_decay_normalized", "subpop")
thetas_shown <- c(1, seq(2, 10, by = 2), seq(25, 100, by = 25))
generations_shown <- seq(0, 50)
migration_shown <- 0

plts <- list()

for (m in seq_along(models_shown)) {
  
  data <- summary_data %>%
    filter(model %in% c("random", models_shown[m]),
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown)
  
  plt <- plot_variable_over_time(data, "parent_corr", "theta", titles,
                                 values, labels) +
    guides(color = guide_legend(byrow = TRUE,
                                override.aes = list(linewidth = lh_linewidth)),
           linetype = "none") +
    theme_manuscript()
  
  lh <- get_legend(plt)
  plt <- plt +
    theme(legend.position = "none")
  
  plts[[m]] <- assemble_panel(plt, title = labels$model[[models_shown[m]]])
  
}

fig <- assemble_figure(plts, lh = lh)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 3.5)

###################################################################

# Figure S8 / Figure S9

models_shown <- c("exponential_decay", "normal_distribution")
thetas_shown <- c(seq(2, 10, by = 2), seq(25, 100, by = 25))
generations_shown <- 20
migration_shown <- 0
seeds_shown <- 2547432215832

for (m in seq_along(models_shown)) {
  fig <- plot_ancestry_hist_by_theta(models_shown[m], thetas_shown,
                                     generations_shown, migration_shown,
                                     seeds_shown)
  counter <- save_figure(figure_directory, fig, counter, 
                         supplementary = TRUE, height = 5.5)
}

#############################################################################

# Figure S10

models_shown <- "subpop"
thetas_shown <- 5
generations_shown <- seq(6, 20, by = 2)
migration_shown <- 0
seeds_shown <- 2547432240285

plts <- list()

for (g in seq_along(generations_shown)) {
  
  sim_shown <- quo(model %in% models_shown &
                   theta %in% thetas_shown &
                   generation %in% generations_shown[g] &
                   migration %in% migration_shown &
                   seed %in% seeds_shown)
  
  rho <- summary_data %>%
    filter(!!sim_shown) %>%
    pull(parent_corr) %>%
    round(digits = 2)

  plt <- binned_ancestry0 %>%
    filter(!!sim_shown, subpop %in% c(0, 1)) %>%
    plot_ancestry_histogram(., "subpop", titles, values$subpop, labels$subpop,
                            alpha = 0.3) +
    labs(title = label_generations(generations_shown[g]),
         subtitle = paste0("r = ", rho)) +
    theme_manuscript() +
    theme(legend.key.height = unit(0.8, "lines"),
          legend.key.width = unit(0.8, "lines"))
  
  y <- layer_scales(plt)$y$range$range[2]
  y <- ceiling(y / 0.1) * 0.1
  ticks <- seq(0, y, length.out = 3)
  
  lh <- get_legend(plt)
  
  plts[[g]] <- plt +
    theme(legend.position = "none")
  
}

panel <- assemble_panel(plts, nrow = 2)
fig <- assemble_figure(list(panel), lh = lh, show_tags = FALSE)
counter <- save_figure(figure_directory, fig, counter, 
                       supplementary = TRUE, height = 3.5)

############################################################################

# Figure S11

models_shown <- c("exponential_decay_normalized", "subpop")
thetas_shown <- seq(4, 10, by = 2)
generations_shown <- c(4, 10, 16, 20)
migration_shown <- 0
seeds_shown <- 2547432215832

panels <- list()

plt_count <- 0

for (p in seq_along(generations_shown)) {
  
  data <- binned_ancestry0 %>%
    filter(model %in% models_shown,
           theta %in% thetas_shown[p],
           generation %in% generations_shown[p],
           migration %in% migration_shown,
           seed %in% seeds_shown, subpop == "all")
  
  y <- max(data$proportion)
  y <- ceiling(y / 0.01) * 0.01
  ticks <- seq(0, y, length.out = 3)
  
  for (m in seq_along(models_shown)) {
    
    theta_title <- labels$theta[as.character(thetas_shown[p])]
    gen_title <- label_generations(generations_shown[p])
    
    plt <- data %>%
      filter(model == models_shown[m]) %>%
      plot_ancestry_histogram(., "model", titles, values$model, labels$model) +
      ggtitle(paste0(theta_title, " (", gen_title, ")")) +
      coord_cartesian(ylim = c(0, y)) +
      scale_y_continuous(breaks = ticks,
                         labels = scales::label_number(accuracy = 0.05)) +
      guides(fill = "none") +
      theme_manuscript() +
      theme(plot.title = element_markdown())
    
    plt_count <- plt_count + 1
    panels[[plt_count]] <- assemble_panel(plt)
  }
  
  }

lh <- binned_ancestry0 %>%
  filter(model %in% models_shown) %>%
  plot_ancestry_histogram(., "model", titles, values$model, labels$model) +
  theme_manuscript() +
  theme(legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(0.8, "lines"))

lh <- get_legend(lh)

fig <- assemble_figure(panels, lh = lh, ncol = 2)

counter <- save_figure(figure_directory, fig, counter, 
                       supplementary = TRUE, height = 6.5)

#############################################################################

# Figure S12

examples_shown <- c(2, 3, 4)
models_shown <- c("exponential_decay_normalized", "subpop")
generations_shown <- 20

panels <- list()

plt_count <- 0

for (e in seq_along(examples_shown)) {
  
  data <- binned_ancestry0 %>%
    filter(model %in% models_shown,
           example %in% examples_shown[e],
           generation %in% generations_shown,
           subpop == "all")
  
  y <- max(data$proportion)
  y <- ceiling(y / 0.01) * 0.01
  ticks <- seq(0, y, length.out = 3)
  
  for (m in seq_along(models_shown)) {
    
    rho <- summary_data %>%
      filter(example == examples_shown[e], model == models_shown[m],
             generation %in% generations_shown) %>%
      pull(parent_corr) %>%
      round(digits = 2)
    
    theta <- summary_data %>%
      filter(example == examples_shown[e], model == models_shown[m],
             generation %in% generations_shown) %>%
      pull(theta)
    
    plt_count <- plt_count + 1
    
    panels[[plt_count]] <- data %>%
      filter(model == models_shown[m]) %>%
      plot_ancestry_histogram(., "model", titles, values$model, labels$model) +
      labs(title = paste0(symbol_theta, " = ", theta),
           subtitle = paste0("r = ", rho)) +
      coord_cartesian(ylim = c(0, y)) +
      scale_y_continuous(breaks = ticks,
                         labels = scales::label_number(accuracy = 0.05)) +
      theme_manuscript() +
      theme(plot.title = element_markdown(),
            legend.position = "none")
  }
  
}

lh <- binned_ancestry0 %>%
  filter(model %in% models_shown) %>%
  plot_ancestry_histogram(., "model", titles, values$model, labels$model) +
  theme_manuscript() +
  theme(legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(0.8, "lines"))

lh <- get_legend(lh)

fig <- assemble_figure(panels, lh = lh, ncol = 2)
counter <- save_figure(figure_directory, fig, counter, 
                       supplementary = TRUE, height = 6.5)

#############################################################################

# Figure S13 / Figure S14 / Figure S15

models_shown <- c("exponential_decay_normalized", "subpop")

fig <- plot_correlation_figure(example_pedigree, geom = "hex",
                               examples = c(2, 4), models = models_shown,
                               n = c(10000, 10000))

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

fig <- plot_correlation_figure(example_pedigree, geom = "point",
                               examples = c(3, 3), models = models_shown,
                               n = c(10000, 100))

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

fig <- plot_correlation_figure(example_pedigree, geom = "hex",
                               examples = c(2, 4), models = models_shown,
                               n = c(100, 100))

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

#############################################################################

# Figure S16

fig <- summary_data %>%
  filter(migration == 0, theta %in% seq(2, 10, 2)) %>%
  ggplot(aes(x = parent_corr, y = ancestry0_var, color = model)) +
  geom_point(size = 0.3) +
  labs(x = titles['parent_corr'], y = titles['ancestry0_var']) +
  guides(color = "none") +
  facet_grid(theta ~ model,
             labeller = labeller(theta = function(x) { 
               paste(symbol_theta, " = ", x) },
               model = labels$short_model )) +
  scale_color_manual(values = values$model, labels = labels$model) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.1)) +
  theme_manuscript() +
  theme(panel.border = element_rect(fill = NA),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_markdown())

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

############################################################################

# Figure S17

models_shown <- c("exponential_decay", "exponential_decay",
                  "normal_distribution", "normal_distribution")
thetas_shown <- 10
generations_shown <- c(15, 20, 10, 12)
migration_shown <- 0
seeds_shown <- 2547432225021

panels <- list()

for (p in seq_along(models_shown)) {
  
  sim_shown <- quo(model == !!models_shown[p] &
                     theta == !!thetas_shown &
                     generation == !!generations_shown[p] &
                     migration == !!migration_shown &
                     seed == !!seeds_shown)
  
  observed <- summary_data %>%
    filter(!!sim_shown) %>%
    pull(parent_corr)
  
  permuted <- shuffled_rho %>%
    filter(!!sim_shown)
  
  plt <- plot_shuffled_rho(observed, permuted) +
    theme_manuscript() +
    guides(color = guide_legend(override.aes = list(linewidth = lh_linewidth)))
  
  y <- max(layer_scales(plt)$y$range$range[2], 0.1)

  lh <- get_legend(plt)
  
  plt <- plt +
    coord_cartesian(ylim = c(0, y)) +
    scale_y_continuous(breaks = seq(0, y, by = 0.05),
                       labels = scales::label_number(accuracy = 0.01)) +
    theme(legend.position = "none")
  
  panels[[p]] <- assemble_panel(plt, title = labels$model[models_shown[p]])
}

fig <- assemble_figure(panels, lh = lh, ncol = 2)

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)


##############################################################################

# Figure S18

examples_shown <- c(1, 3, 3)
models_shown <- c("random", "exponential_decay_normalized", "subpop")

x <- example_pedigree %>%
  filter(example %in% examples_shown) %>%
  summarize(lim1 = min(across(parent1_ancestry0:parent2_ancestry0)),
            lim2 = max(across(parent1_ancestry0:parent2_ancestry0)))

x <- c(floor(x$lim1 / 0.05) * 0.05, ceiling(x$lim2 / 0.05) * 0.05)
xticks <- seq(x[1], x[2]-0.05, length.out = 4)

panels <- list()

for (e in seq_along(examples_shown)) {
  
  plt <- example_pedigree %>%
    filter(example == examples_shown[e], model == models_shown[e]) %>%
    plot_mating_pairs_hurricane(.) +
    theme_manuscript() +
    labs(x = titles$ancestry0, y = "Individual") +
    coord_cartesian(xlim = c(x[1], x[2])) +
    scale_x_continuous(breaks = xticks,
                       labels = scales::label_number(accuracy = 0.1)) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  panels[[e]] <- assemble_panel(plt,
                                title = labels$short_model[models_shown[e]])
  
}

fig <- assemble_figure(panels, ncol = 3)

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 2.5)

############################################################################

# Figure S19

models_shown <- c("exponential_decay", "exponential_decay",
                  "normal_distribution", "normal_distribution")
thetas_shown <- 10
generations_shown <- c(15, 20, 10, 12)
migration_shown <- 0
seeds_shown <- 2547432225021

panels <- list()

for (p in seq_along(models_shown)) {
  
  sim_shown <- quo(model == !!models_shown[p] &
                     theta == !!thetas_shown &
                     generation == !!generations_shown[p] &
                     migration == !!migration_shown &
                     seed == !!seeds_shown)
  
  observed <- summary_data %>%
    filter(!!sim_shown) %>%
    pull(delta_obs)
  
  permuted <- shuffled_delta %>%
    filter(!!sim_shown)
  
  plt <- plot_shuffled_delta(observed, permuted) +
    theme_manuscript2() +
    guides(color = guide_legend(override.aes = list(linewidth = lh_linewidth)))
  
  y <- max(layer_scales(plt)$y$range$range[2], 0.1)
  
  lh <- get_legend(plt)
  
  plt <- plt +
    coord_cartesian(ylim = c(0, y)) +
    scale_y_continuous(breaks = seq(0, y, by = 0.05),
                       labels = scales::label_number(accuracy = 0.01)) +
    theme(legend.position = "none")
  
  panels[[p]] <- assemble_panel(plt, title = labels$model[models_shown[p]])
}

fig <- assemble_figure(panels, lh = lh, ncol = 2)

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

############################################################################

# Figure S20

models_shown <- models
thetas_shown <- seq(1, 10)
generations_shown <- seq(3, 18, by = 3)
migration_shown <- 0

data <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown)

y <- max(data$parent_corr)
y <- ceiling(y / 0.1) * 0.1
ticks <- seq(0, y, length.out = 4)

plts <- list()

for (g in seq_along(generations_shown)) {
  plt <- data %>%
    filter(generation == generations_shown[g]) %>%
    plot_correlation_vs_B(., titles, values, labels) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                override.aes = list(size = lh_point_size))) +
    coord_cartesian(ylim = c(0, y)) +
    scale_y_continuous(breaks = ticks,
                       labels = scales::label_number(accuracy = 0.1)) +
    theme_manuscript()
  
  lh <- get_legend(plt)
  
  plts[[g]] <- plt +
    theme(legend.position = "none")
  
}

panel <- assemble_panel(plts, nrow = 2)
fig <- assemble_figure(list(panel), lh = lh, show_tags = FALSE)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5)

##############################################################################

# Figure S21

models_shown <- c("exponential_decay_normalized", "subpop")
thetas_shown <- seq(2, 10)
generations_shown <- seq(3, 18, by = 3)
migration_shown <- 0

data <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown)

y <- max(data$delta_obs)
y <- ceiling(y / 0.1) * 0.1
ticks <- seq(0, y, length.out = 3)

panels <- list()

for (g in seq_along(generations_shown)) {
  
  plt <- data %>%
    filter(generation == generations_shown[g]) %>%
    plot_correlation_vs_delta(., titles, values, labels) +
    guides(color = guide_legend(override.aes = list(size = lh_point_size))) +
    coord_cartesian(ylim = c(0, y)) +
    scale_y_continuous(breaks = ticks,
                       labels = scales::label_number(accuracy = 0.1))
  
  lh <- get_legend(plt)
  
  plt <- plt +
    theme(legend.position = "none")
  
  panels[[g]] <- assemble_panel(plt)
  
}

fig <- assemble_figure(panels, lh = lh, show_tags = FALSE)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 4.5)

##############################################################################

# Figure S22

models_shown <- c("exponential_decay", "exponential_decay_normalized",
                  "normal_distribution", "subpop")
thetas_shown <- seq(2, 10, by = 1)
generations_shown <- 20
migration_shown <- 0

panels <- list()

for (m in seq_along(models_shown)) {
  
  plt <- summary_data %>%
    filter(model == models_shown[m],
           theta %in% thetas_shown,
           generation %in% generations_shown,
           migration %in% migration_shown) %>%
    ggplot(aes(x = theta, y = tracts_q50, color = model)) +
    geom_point() +
    theme_manuscript() +
    labs(x = titles['theta'], y = titles['tract_median']) +
    scale_color_manual(values = values$model, labels = labels$model) +
    guides(color = "none") +
    theme(plot.title = element_text(face = "bold"),
          axis.title.x = element_markdown()) +
    coord_cartesian(ylim = c(6, 9))
  
  panels[[m]] <- assemble_panel(plt, title = labels$model[models_shown[m]])
}

fig <- assemble_figure(panels)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5.5)

##############################################################################

# Figure S23

examples_shown <- c(2, 3, 4)
generations_shown <- seq(2, 50)

panels <- list()

for (e in seq_along(examples_shown)) {
  
    rho <- summary_data %>%
      filter(example == examples_shown[e] & generation == 20) %>%
      pull(parent_corr) %>%
      round(digits = 2) %>%
      unique() %>%
      sprintf(fmt = "%0.02f", .)
    
    data <- summary_data %>%
      filter(example == examples_shown[e],
             generation %in% generations_shown)
    
    y <- max(data$timing_discrepancy)
    y <- ceiling(y)
    if (y > 2) { ticks <- seq(0, y, by = 2) } else { ticks <- seq(0, y) }
    
    plt <- data %>%
      plot_variable_over_time(., "timing_discrepancy", "model", titles, values,
                              labels) +
      labs(title = paste0("r = ", rho)) +
      coord_cartesian(ylim = c(0, y)) +
      scale_y_continuous(breaks = ticks,
                         labels = scales::label_number(accuracy = 1)) +
      guides(color = guide_legend(override.aes =
                                    list(linewidth = lh_linewidth)),
             linetype = "none") +
      theme_manuscript()
    
    lh <- get_legend(plt)
    
    plt <- plt +
      theme(legend.position = "none")
    
    panels[[e]] <- assemble_panel(plt)
    
  }

fig <- assemble_figure(panels, lh = lh, ncol = 3) 
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 2.5)

##############################################################################

# Figure S24

models_shown <- models
thetas_shown <- seq(1, 10)
generations_shown <- c(5, 20)
migration_shown <- 0

data <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_rates)

x <- max(data$parent_corr)
x <- ceiling(x / 0.1) * 0.1
xticks <- seq(0, x, length.out = 5)[1:4]

y <- max(data$timing_discrepancy)
y <- ceiling(y)
if (y > 2) { yticks <- seq(0, y, by = 2) } else { yticks <- seq(0, y) }

panels <- list()

for (g in seq_along(generations_shown)) {
  
  plt <- data %>%
    filter(generation == generations_shown[g]) %>%
    plot_discrepancy_in_timing(., values, labels) +
    theme_manuscript() +
    labs(x = paste(titles['parent_corr'], "at",
                   label_generations(generations_shown[g], short = TRUE)),
         y = paste(titles['timing_discrepancy'], "at Gen. 20")) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                override.aes = list(size = lh_point_size))) +
    coord_cartesian(xlim = c(-0.05, x), ylim = c(-0.5, y)) +
    scale_x_continuous(breaks = xticks,
                       labels = scales::label_number(accuracy = 0.1)) +
    scale_y_continuous(breaks = yticks,
                       labels = scales::label_number(accuracy = 1))
  
  lh <- get_legend(plt)
  panels[[g]] <- plt +
    theme(legend.position = "none")
}

fig <- assemble_figure(panels, ncol = 2, lh = lh)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 3.5)

##############################################################################

# Figure S25

models_shown <- c("exponential_decay", "normal_distribution")
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
    guides(color = guide_legend(override.aes =
                                  list(linewidth = lh_linewidth)),
           linetype = guide_legend(override.aes =
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
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5.5)
##############################################################################

# Figure S26

models_shown <- models
thetas_shown <- thetas
generations_shown <- seq(0, 20)
migration_shown <- c(0, 0.01)

data <- summary_data %>%
  filter(model %in% models_shown,
         theta %in% thetas_shown,
         generation %in% generations_shown,
         migration %in% migration_shown) %>%
  group_by(generation, model, theta, migration) %>%
  summarize(rho = mean(parent_corr)) %>%
  pivot_wider(names_from = "migration", names_prefix = "migration_",
              values_from = "rho") %>%
  mutate(delta = migration_0 - migration_0.01)

annotations <- tibble(x = 1, y = c(0.3, -0.3),
                      lab = c("Correlation greater under single-pulse admixture",
                              "Correlation greater under continuous-migration"))

fig <- ggplot() +
  geom_path(data = data,
            mapping = aes(x = generation, y = delta, color = model,
                          group = interaction(model, theta))) +
  geom_text(data = annotations,
            mapping = aes(x = x, y = y, label = lab),
            family = "Avenir Next", hjust = 0) +
  geom_path(aes(x = seq(0.5, 20.5), y = 0), linetype = "dotted",
            linewidth = 1) +
  scale_color_manual(values = values$model, labels = labels$model) +
  labs(x = titles['generation'],
       y = "Change in correlation between scenarios") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = lh_linewidth))) +
  coord_cartesian(ylim = c(-0.35, 0.35)) +
  theme_manuscript()

counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5.5)
##############################################################################

# Figure S27

models_shown <-  c("exponential_decay_normalized", "subpop",
                   "exponential_decay", "normal_distribution")
thetas_shown <- c(1, 2, 4, 10)
generations_shown <- seq(0, 50)
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
    guides(color = guide_legend(override.aes =
                                  list(linewidth = lh_linewidth)),
           linetype = guide_legend(override.aes =
                                     list(linewidth = lh_linewidth))) +
    theme_manuscript() +
    theme(legend.box = "vertical")
  
  lh <- get_legend(plt)
  
  plt <- plt + theme(legend.position = "none")
  
  panels[[m]] <- assemble_panel(plt) +
    add_title(labels$model[models_shown[m]])
  
}

fig <- assemble_figure(panels, lh = lh, ncol = 2)
counter <- save_figure(figure_directory, fig, counter, supplementary = TRUE,
                       height = 5.5)
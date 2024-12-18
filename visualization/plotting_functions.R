library(ggplot2)

plot_ancestry_histogram <- function(data, group, titles, colors, labels,
                                    alpha = 1) {
  
  plt <- ggplot(data = data,
                mapping = aes(x = ancestry0, y = proportion,
                              fill = .data[[group]])) +
    geom_col(width = 0.02, position = "identity", alpha = alpha) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                       labels = c(0, 0.25, 0.50, 0.75, 1)) +
    scale_fill_manual(values = colors, labels = labels) +
    labs(x = titles["ancestry0"],
         y = titles["proportion"]) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
  
  y <- layer_scales(plt)$y$range$range[2]
  y <- ceiling(y / 0.1) * 0.1
  ticks <- seq(0, y, length.out = 3)
  
  plt <- plt +
    coord_cartesian(default = TRUE, ylim = c(0, y)) +
    scale_y_continuous(breaks = ticks,
                       labels = scales::label_number(accuracy = 0.01))
  return(plt)
  
}

plot_psi_over_time <- function(data, generations, colors, labels) {
  
  plt <- ggplot(data = data,
                mapping = aes(x = ancestry0, y = psi_ratio,
                              color = as.factor(generation))) +
    geom_point(size = 0.01) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                       labels = c(0, 0.25, 0.50, 0.75, 1)) +
    scale_color_manual(values = colors,
                       labels = label_generations(generations)) +
    labs(x = bquote(x[i]),
         y = paste0("max(", symbol_psi_i, ") / min(", symbol_psi_i, ")")) +
    guides(color = guide_legend(byrow = TRUE,
                                override.aes = list(size = 1.2)))
  
  return(plt)
}

plot_offspring_vs_ancestry <- function(data, control) {
  
  plt <- ggplot() +
    geom_ribbon(data = control,
                aes(x = ancestry0, ymin = ymin, ymax = ymax),
                alpha = 0.2) +
    geom_point(data = data,
               aes(x = ancestry0, y = proportion, color = parent)) +
    scale_color_manual(values = c("#e41a1c", "#377eb8"),
                       labels = c("Parent 1", "Parent 2")) +
    labs(x = titles["ancestry0"],
         y = titles["proportion"],
         title = labels$model[unique(as.character(data$model))])
  
  return(plt)
}


plot_variable_over_time <- function(data, variable, key, titles, colors, labels) {
  
  plt <- ggplot(data = data,
                mapping = aes(x = generation, y = .data[[variable]],
                              color = .data[[key]], linetype = migration,
                              group = simulation)) +
    geom_path() +
    scale_color_manual(values = colors[[key]], labels = labels[[key]]) +
    scale_linetype_manual(values = colors$migration,
                          labels = labels$migration) +
    labs(x = titles['generation'], y = titles[[variable]]) +
    guides(color = guide_legend(byrow = TRUE))
  
  return(plt)
  
}

plot_correlation_vs_delta <- function(data, titles, values, labels) {
  
  g <- unique(data$generation)
  
  plt <- ggplot(data = data,
                mapping = aes(x = parent_corr, y = delta_obs, color = model)) +
    geom_point() +
    labs(x = titles['parent_corr'], y = titles['obs_delta']$obs_delta,
         title = label_generations(g)) +
    scale_color_manual(values = values$model, labels = labels$model) +
    theme_manuscript()
  
  return(plt)
  
}

plot_parent_correlation <- function(data, geom, ...) {
  
  lims <- c(0.3, 0.7)
  
  rho <- sprintf("%0.2f", cor(data$parent1_ancestry0, data$parent2_ancestry0))
  n <- format(nrow(data), big.mark = ",")
  n <- tibble(x = 0.32, y = 0.3, lab = paste0("n = ", n, " individuals"))
  
  plt <- ggplot(data = data,
                mapping = aes(x = parent1_ancestry0, y = parent2_ancestry0)) +
    labs(x = titles$parent1_ancestry0, y = titles$parent2_ancestry0,
         subtitle = bquote(.(titles$parent_corr)~"="~.(rho))) +
    coord_cartesian(xlim = lims, ylim = lims) +
    scale_x_continuous(breaks = seq(lims[1], lims[2], 0.1)) +
    scale_y_continuous(breaks = seq(lims[1], lims[2], 0.1))
  
  if (geom == "point") {
    plt <- plt +
      geom_point(...)
  } else if (geom == "hex") {
    plt <- plt +
      geom_hex(mapping = aes(fill = after_stat(ndensity)),
               binwidth = 0.025) +
      scico::scale_fill_scico(256, palette = "davos", direction = -1) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                                   nbin = 8, frame.colour = "#000000",
                                   frame.linewidth = 0.3,
                                   theme = theme(legend.ticks = element_blank())))
  }
  
  plt <- plt +
    geom_abline(slope = 1, intercept = 0, color = "#BDBDBD", linewidth = 0.3) +
    geom_text(data = n, aes(x = x, y = y, label = lab),
              family = "Avenir Next", size = 8/.pt, hjust = "left")
  
  return(plt)
  
}


plot_correlation_figure <- function(pedigree, geom, examples, models, n) {
  
  panels <- list()
  
  counter <- 0
  
  for (p in seq_along(examples)) {
    
    for (m in seq_along(models)) {
      
      counter <- counter + 1
      
      set.seed(n[p])
      data <- pedigree %>%
        filter(example == examples[p], model == models[m]) %>%
        sample_n(n[p], replace = FALSE)
      
      plt <- plot_parent_correlation(data, geom) +
        theme_manuscript2() +
        theme(legend.position = "inside",
              legend.position.inside = c(0.97, 0.35))
      
      if (counter > 1) { plt <- plt + theme(legend.position = "none") }
      
      title <- labels$model[unique(as.character(data$model))]
      panels[[counter]] <- assemble_panel(plt) + add_title(title)
      
    }
  }
  
  figure <- assemble_figure(panels, nrow = 2)
  return(figure)
  
}

plot_ancestry_psi <- function(data, titles, values, labels) {
  
  plt <- ggplot(data = data,
                mapping = aes(x = xi, y = psi, color = theta,
                              linetype = generation)) +
    geom_path() +
    scale_color_manual(values = values$theta, labels = labels$theta) +
    scale_linetype_manual(values = values$generation,
                          labels = label_generations(seq(3))) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                       labels = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(x = titles$xi, y = titles$psi,
         # title = labels$short_model[as.character(unique(data$model))],
         subtitle = bquote(x[i] == .(x1_shown[x])))
  
  if (unique(data$model) == "exponential_decay_normalized") {
    plt <- plt +
      guides(color = "none")
  } else {
    plt <- plt +
      guides(linetype = "none")
  }
  
  return(plt)
}

plot_social_psi <- function(data, titles, values, labels) {
  plt <- ggplot(data = data,
                mapping = aes(x = social_group, y = psi, fill = theta)) +
    geom_col(position = "dodge") +
    scale_x_discrete(labels = c("same" = bquote(s[i]==s[j]),
                                "different" = bquote(s[i]!=s[j]))) +
    scale_fill_manual(values = values$theta, labels = labels$theta) +
    labs(x = titles['social_group'], y = titles$psi)
  
  return(plt)
}

plot_correlation_vs_B <- function(data, titles, values, labels) {
  
  g <- unique(data$generation)
  
  plt <- ggplot(data = data,
                mapping = aes(x = B, y = parent_corr,
                              color = model)) +
    geom_point() +
    scale_color_manual(values = values$model,
                       labels = labels$model) +
    labs(x = titles["expressed_mating_bias"],
         y = titles["parent_corr"],
         title = label_generations(g))
  
  return(plt)
}

plot_ancestry_histogram_over_time <- function(simulation_id, generations_shown) {
  
  plts <- list()
  for (g in generations_shown) {
    
    v <- summary_data %>%
      filter(simulation == simulation_id, generation == g) %>%
      pull(ancestry0_var) %>%
      sprintf(fmt = "%0.04f", .)
    
    data <- binned_ancestry0 %>%
      filter(simulation == simulation_id, generation == g, subpop == "all")
    
    plt <-  data %>%
      plot_ancestry_histogram(., "model", titles, values$model, labels$model) +
      guides(fill = "none") +
      labs(title = paste0("t = ", g),
           subtitle = paste0(symbol_sigma, " = ", v)) +
      theme_manuscript() +
      theme(plot.subtitle = element_markdown(size = 6))
    
    plts[[g]] <- plt
    
  }
  
  fig <- cowplot::plot_grid(plotlist = plts, ncol = 4)
  
  title <- labels$model[unique(data$model)]
  fig <- fig + add_title(title)
  
  return(fig)
}

plot_estimated_vs_true_time <- function(example) {
  
  rho <- summary_data %>%
    filter(example == !!example, generation == 20) %>%
    pull(parent_corr) %>%
    round(digits = 2) %>%
    unique() %>%
    sprintf(fmt = "%0.02f", .)
  
  plt <- summary_data %>%
    filter(example %in% c(1, !!example), generation > 1, generation <= 20) %>%
    plot_variable_over_time(., "timing_0", "model", titles,
                            values, labels) +
    labs(x = "True t", title = bquote(.(titles$parent_corr)~"="~.(rho))) +
    theme_manuscript() +
    guides(linetype = "none")
  
  return(plt)
  
}



plot_ancestry_hist_by_theta <- function(model, thetas, generation, migration,
                                        initial, seed) {
  
  plts <- list()
  
  for (t in seq_along(thetas)) {
    
    sim_shown <- quo(model == !!model & theta == thetas[t] &
                       seed == !!seed & migration == !!migration &
                       initial == !!initial &
                       generation == !!generation)
    
    data <- binned_ancestry0 %>%
      filter(!!sim_shown, subpop == "all")
    
    title <- labels$theta[as.character(unique(data$theta))]
    
    rho <- summary_data %>%
      filter(!!sim_shown) %>%
      pull(parent_corr) %>%
      round(., digits = 2)
    
    rho <- sprintf(fmt = "%0.02f", rho + 0)
    
    plts[[t]] <- plot_ancestry_histogram(data, "model", titles, values$model,
                                         labels$model) +
      labs(title = title, subtitle = bquote(.(titles$parent_corr)~"="~ .(rho))) +
      guides(fill = "none") +
      theme_manuscript() +
      theme(plot.title = element_markdown())
    
  }
  
  fig <- assemble_panel(plts, title = labels$model[unique(data$model)])
  return(fig)
  
}

plot_discrepancy_in_timing <- function(data, generation, values, labels) {
  
  rho <- cor(data[[as.name(paste0('parent_corr_', generation))]],
             data$timing_discrepancy_20)
  rho <- sprintf(fmt = "%0.02f", rho)
  
  plt <- data %>%
    ggplot(aes(x =.data[[paste0('parent_corr_', generation)]],
               y = timing_discrepancy_20, color = model)) +
    geom_point() +
    scale_color_manual(values = values$model, labels = labels$model) +
    labs(title = paste0("r = ", rho))
  
  return(plt)
}

plot_mating_pairs_hurricane <- function(data, ...) {
  
  data <- data %>%
    arrange(delta) %>%
    mutate(n = seq(nrow(.)))
  
  plt <- ggplot(data = data,
                aes(x = parent1_ancestry0, xend = parent2_ancestry0, y = n)) +
    geom_segment(...)
  
  return(plt)
  
}

plot_permutation_test <- function(observed, permuted, var, pval, xlab, values,
                                  labels) {
  
  generation <- unique(permuted$generation)
  
  summary <- tibble(observed = observed,
                    permuted = mean(permuted[[var]])) %>%
    pivot_longer(cols = c("permuted", "observed"), names_to = "x",
                 values_to = "mu")
  
  plt <- ggplot() +
    geom_histogram(data = permuted,
                   mapping = aes(x = .data[[var]],
                                 y = after_stat(count) /
                                   sum(after_stat(count))),
                   fill = "#CCCCCC", alpha = 0.5, bins = 30) +
    geom_vline(data = summary, mapping = aes(xintercept = mu, color = x),
               key_glyph = "path") +
    labs(x = xlab, y = "Proportion of simulations",
         subtitle = paste0(label_generations(generation),
                           " (", pval, ")")) +
    scale_color_manual(values = values, labels = labels) 
  
  return(plt)
}

plot_shuffled_delta <- function(observed, permuted) {
  
  pval <- (sum(permuted$delta_perm <= observed)+1) / (nrow(permuted)+1)
  
  pval <- paste0("p = ", sprintf(fmt = "%0.03f", pval))
  s <- case_when(pval < 1/1000 ~ "<", .default = "=")
  pval <- paste("p", s, sprintf(fmt = "%0.03f", pval))
  
  values <- c("observed" = "#2C7FB8", "permuted" = "#636363")
  labels <- c("observed" = bquote(bar(Delta)[obs]),
              "permuted" = bquote(bar(Delta)[perm]))
  xlab <- bquote(bar(Delta[x]))
  
  plt <- plot_permutation_test(observed, permuted, "delta_perm", pval,
                               xlab, values, labels)
  return(plt)
}

plot_shuffled_rho <- function(observed, permuted, titles) {
  
  pval <- (sum(permuted$rho_perm >= observed)+1) / (nrow(permuted)+1)
  s <- case_when(pval < 1/1000 ~ "<", .default = "=")
  pval <- paste("p", s, sprintf(fmt = "%0.03f", pval))
  
  values <- c("observed" = "#2C7FB8", "permuted" = "#636363")
  labels <- c("observed" = bquote(Observed~.(titles$parent_corr)),
              "permuted" = bquote(Mean~permuted~.(titles$parent_corr)))
  
  xlab <- titles$parent_corr
  
  plt <- plot_permutation_test(observed, permuted, "rho_perm", pval,
                               xlab, values, labels)
  
  
  return(plt)
}
